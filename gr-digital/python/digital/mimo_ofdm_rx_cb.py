#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2018 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr
import digital_swig as digital
from gnuradio.digital.mimo_decoder_cc import mimo_decoder_cc

try:
    # This will work when feature #505 is added.
    from gnuradio import fft
    from gnuradio import blocks
    from gnuradio import analog
except ImportError:
    # Until then this will work.
    import fft_swig as fft
    import blocks_swig as blocks
    import analog_swig as analog

import numpy as np
from scipy.linalg import hadamard

''' 
Default values for Receiver.
'''
# Number of receiving antennas.
_def_n = 2
_def_mimo_technique = 'vblast'
_def_fft_len = 64
_def_cp_len = 16
_seq_seed = 42
# Modulation order of header and payload symbols
_def_bps_header = 1
_def_bps_payload = 1
# Tag keys.
_def_start_key = "start"
_def_csi_key = "csi"
_def_frame_length_tag_key = "frame_length"
_def_packet_length_tag_key = "packet_length"
_def_packet_num_tag_key = "packet_num"
# OFDM carrier allocation.
_def_pilot_carriers_mimo=[range(-26, 0, 3)+range(2, 27, 3), ]
_def_occupied_carriers_mimo = [[x for x in range(-24, 25, 1) if x not in _def_pilot_carriers_mimo[0] + [0]], ]

def _get_active_carriers(fft_len, occupied_carriers, pilot_carriers):
    """ Returns a list of all carriers that at some point carry data or pilots. """
    active_carriers = list()
    for carrier in list(occupied_carriers[0]) + list(pilot_carriers[0]):
        if carrier < 0:
            carrier += fft_len
        active_carriers.append(carrier)
    return active_carriers

def _make_sync_word1(fft_len, occupied_carriers, pilot_carriers):
    """ Creates a random sync sequence for fine frequency offset and timing
    estimation. This is the first of typically two sync preamble symbols
    for the Schmidl & Cox sync algorithm.
    The relevant feature of this symbols is that every second sub-carrier
    is zero. In the time domain, this results in two identical halves of
    the OFDM symbols.
    Symbols are always BPSK symbols. Carriers are scaled by sqrt(2) to keep
    total energy constant.
    Carrier 0 (DC carrier) is always zero. If used, carrier 1 is non-zero.
    This means the sync algorithm has to check on odd carriers!
    """
    active_carriers = _get_active_carriers(fft_len, occupied_carriers, pilot_carriers)
    np.random.seed(_seq_seed)
    bpsk = {0: np.sqrt(2), 1: -np.sqrt(2)}
    sw1 = [bpsk[np.random.randint(2)]  if x in active_carriers and x % 2 else 0 for x in range(fft_len)]
    return np.fft.fftshift(sw1)


def _make_sync_word2(fft_len, occupied_carriers, pilot_carriers):
    """ Creates a random sync sequence for coarse frequency offset and channel
    estimation. This is the second of typically two sync preamble symbols
    for the Schmidl & Cox sync algorithm.
    Symbols are always BPSK symbols.
    """
    active_carriers = _get_active_carriers(fft_len, occupied_carriers, pilot_carriers)
    np.random.seed(_seq_seed)
    bpsk = {0: 1, 1: -1}
    sw2 = [bpsk[np.random.randint(2)] if x in active_carriers else 0 for x in range(fft_len)]
    sw2[0] = 0j
    return np.fft.fftshift(sw2)

def _get_constellation(bps):
    """ Returns a modulator block for a given number of bits per symbol """
    constellation = {
            1: digital.constellation_bpsk(),
            2: digital.constellation_qpsk(),
            3: digital.constellation_8psk()
    }
    try:
        return constellation[bps]
    except KeyError:
        print 'Modulation not supported.'
        exit(1)

class mimo_ofdm_rx_cb(gr.hier_block2):
    """
    Hierarchical block for MIMO-OFDM demodulation.

    The inputs are N complex baseband signals (e.g. from UHD sources).
    The detected packets are output as a stream of packed bits on the output.

    Args:
    n: Number of receiving antennas.
    mimo_technique: Key to one of the implemented MIMO techniques. Choose between:
    |               V-BLAST (key = 'vblast')
    |               Alamouti (key = 'alamouti')
    |               Differential Alamouti (key = 'diff_stbc')
    fft_len: The length of FFT (integer).
    cp_len:  The length of cyclic prefix in total samples (integer).
    frame_length_tag_key: Used internally to tag the length of the OFDM frame.
    packet_length_tag_key: The name of the tag giving packet length at the input.
    occupied_carriers: A vector of vectors describing which OFDM carriers are occupied.
    pilot_carriers: A vector describing which OFDM carriers are occupied with pilot symbols.
    pilot_symbols: The pilot symbols.
    bps_header: Bits per symbol (header).
    bps_payload: Bits per symbol (payload).
    sync_word1: The first sync preamble symbol. This has to be with
    |           zeros on alternating carriers. Used for fine and
    |           coarse frequency offset and timing estimation.
    sync_word2: The second sync preamble symbol. This has to be filled
    |           entirely. Also used for coarse frequency offset and
    |           channel estimation.
    """
    def __init__(self,
                 n=_def_n,
                 mimo_technique=_def_mimo_technique,
                 fft_len=_def_fft_len,
                 cp_len=_def_cp_len,
                 start_key=_def_start_key,
                 csi_key = _def_csi_key,
                 frame_length_tag_key=_def_frame_length_tag_key,
                 packet_length_tag_key=_def_packet_length_tag_key,
                 packet_num_tag_key=_def_packet_num_tag_key,
                 occupied_carriers=_def_occupied_carriers_mimo,
                 pilot_carriers=_def_pilot_carriers_mimo,
                 pilot_symbols=None,
                 bps_header=_def_bps_header,
                 bps_payload=_def_bps_payload,
                 sync_word1=None,
                 sync_word2=None,
                 scramble_bits=False):
        gr.hier_block2.__init__(self,
            "mimo_ofdm_rx_cb",
            gr.io_signature(n, n, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature2(2, 2, gr.sizeof_char, gr.sizeof_gr_complex))  # Output signature

        """
        Parameter initalization
        """
        self.n = n
        self.mimo_technique = mimo_technique
        self.fft_len = fft_len
        self.cp_len = cp_len
        self.start_key = start_key
        self.csi_key = csi_key
        self.frame_length_tag_key = frame_length_tag_key
        self.packet_length_tag_key = packet_length_tag_key
        self.occupied_carriers = occupied_carriers
        self.pilot_carriers = pilot_carriers,
        self.bps_header = bps_header
        self.bps_payload = bps_payload

        # Change SISO/MIMO specific default parameters.
        if occupied_carriers is None:
            self.occupied_carriers = _def_occupied_carriers_mimo
        if pilot_carriers is None:
            self.pilot_carriers = _def_pilot_carriers_mimo

        if pilot_symbols is None:
            # Generate Hadamard matrix as orthogonal pilot sequences.
            self.pilot_symbols = hadamard(_def_n)
        else:
            self.pilot_symbols = pilot_symbols



        if sync_word1 is None:
            self.sync_word1 = _make_sync_word1(fft_len, occupied_carriers, pilot_carriers)
        else:
            if len(sync_word1) != self.fft_len:
                raise ValueError("Length of sync sequence(s) must be FFT length.")
            self.sync_word1 = sync_word1
        self.sync_word2 = ()
        if sync_word2 is None:
            self.sync_word2 = _make_sync_word2(fft_len, occupied_carriers, pilot_carriers)
        elif len(sync_word2):
            if len(sync_word2) != fft_len:
                raise ValueError("Length of sync sequence(s) must be FFT length.")
            self.sync_word2 = sync_word2
        if scramble_bits:
            self.scramble_seed = 0x7f
        else:
            self.scramble_seed = 0x00 # We deactivate the scrambler by init'ing it with zeros

        """
        Synchronization (timing sync and fractional carrier frequency sync) 
        """
        add = blocks.add_cc()
        sum_sync_detect = digital.ofdm_sync_sc_cfb(fft_len, cp_len)
        mimo_sync = digital.mimo_ofdm_synchronizer_fbcvc(self.n,
                                                         self.fft_len,
                                                         self.cp_len,
                                                         self.sync_word1,
                                                         self.sync_word2,
                                                         self.start_key)
        # Factor for OFDM energy normalization.
        rx_normalize = 1.0 / np.sqrt(self.fft_len)
        symbol_len = fft_len + cp_len
        manual_adjusting_factor = 3 + cp_len
        for i in range(0, self.n):
            # Add up MIMO signals to do the sync on this reference signal. #TODO delays set properly???
            self.connect((self, i), blocks.multiply_const_cc(rx_normalize), (add, i))
            self.connect((self, i), blocks.multiply_const_cc(rx_normalize),
                        blocks.delay(gr.sizeof_gr_complex, symbol_len),
                        (mimo_sync, 3+i))
        self.connect(add, sum_sync_detect)
        self.connect((sum_sync_detect, 0), (mimo_sync, 0))  # Fine frequency offset signal.
        self.connect((sum_sync_detect, 1), blocks.delay(gr.sizeof_char, manual_adjusting_factor), (mimo_sync, 1))  # Trigger signal.
        self.connect(add, blocks.delay(gr.sizeof_gr_complex, symbol_len), (mimo_sync, 2))  # Sum signal.

        """
        OFDM demodulation
        """
        ofdm_demod = []
        for i in range(0, self.n):
            ofdm_demod.append(fft.fft_vcc(self.fft_len, True, (), True))
            self.connect((mimo_sync, i), ofdm_demod[i])
        # self.connect(add_tag1, ofdm_demod[0])
        # self.connect(add_tag2, ofdm_demod[1])

        """
        Carrier frequency correction
        """
        carrier_freq_corrector = digital.ofdm_correct_carrier_freq_offset_vcvc(
            self.n, self.fft_len, self.cp_len, "carrier_freq_offset")
        for i in range(0, self.n):
            self.connect(ofdm_demod[i], (carrier_freq_corrector, i))

        """
        MIMO channel estimation
        """
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(
            n=self.n,
            fft_len=fft_len,
            pilot_symbols=self.pilot_symbols,
            pilot_carriers=pilot_carriers[0],
            occupied_carriers=self.occupied_carriers[0],
            csi_key=self.csi_key,
            start_key=self.start_key)
        for i in range(0, self.n):
            self.connect((carrier_freq_corrector, i), (channel_est, i))

        """
        MIMO decoder
        """
        mimo_decoder = mimo_decoder_cc(
            N=self.n,
            mimo_technique=self.mimo_technique,
            vlen=len(self.occupied_carriers[0]),
            csi_key=self.csi_key)
        for i in range(0, self.n):
            self.connect((channel_est, i), (mimo_decoder, i))

        """
        Header reader/parser
        """
        header_constellation = _get_constellation(bps_header)
        header_formatter = digital.packet_header_ofdm(
            occupied_carriers, 1,
            packet_length_tag_key,
            frame_length_tag_key,
            packet_num_tag_key,
            bps_header, bps_payload)
        header_reader = digital.mimo_ofdm_header_reader_cc(header_constellation.base(),
                                                           header_formatter.formatter(),
                                                           self.start_key)
        self.connect(mimo_decoder, header_reader)

        """
        Payload demodulation + CRC
        """
        payload_constellation = _get_constellation(bps_payload)
        payload_demod = digital.constellation_decoder_cb(payload_constellation.base())
        payload_pack = blocks.repack_bits_bb(bps_payload, 8, self.packet_length_tag_key, True)
        crc = digital.crc32_bb(True, self.packet_length_tag_key)
        self.connect(header_reader, payload_demod, payload_pack, crc, self)
        self.connect(mimo_decoder, (self, 1))
