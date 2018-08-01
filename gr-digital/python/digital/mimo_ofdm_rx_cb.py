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

_def_n = 2
_def_mimo_technique = 'none'
_def_fft_len = 64
_def_cp_len = 16
_def_frame_length_tag_key = "frame_length"
_def_packet_length_tag_key = "packet_length"
_def_packet_num_tag_key = "packet_num"
# Data and pilot carriers are same as in 802.11a
_def_occupied_carriers = (range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22, 27),)
_def_pilot_carriers=((-21, -7, 7, 21,),)
_pilot_sym_scramble_seq = (
        1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1,
        1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1,
        -1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1,
        -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1
)
_def_pilot_symbols= tuple([(x, x, x, -x) for x in _pilot_sym_scramble_seq])
_seq_seed = 42

_walsh_sequences = np.array([
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1,-1, 1,-1, 1,-1, 1,-1],
    [1, 1,-1,-1, 1, 1,-1,-1],
    [1,-1,-1, 1, 1,-1,-1, 1],
    [1, 1, 1, 1,-1,-1,-1,-1],
    [1,-1, 1,-1,-1, 1,-1, 1],
    [1, 1,-1,-1,-1,-1, 1, 1],
    [1,-1,-1, 1,-1, 1, 1,-1]
])


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
    docstring for block mimo_ofdm_rx_cb
    """
    def __init__(self, n=_def_n, mimo_technique=_def_mimo_technique,
                 fft_len=_def_fft_len, cp_len=_def_cp_len,
                 frame_length_tag_key=_def_frame_length_tag_key,
                 packet_length_tag_key=_def_packet_length_tag_key,
                 packet_num_tag_key=_def_packet_num_tag_key,
                 occupied_carriers=_def_occupied_carriers,
                 pilot_carriers=_def_pilot_carriers,
                 pilot_symbols=_walsh_sequences,
                 bps_header=1,
                 bps_payload=1,
                 sync_word1=None,
                 sync_word2=None,
                 debug_log=False,
                 scramble_bits=False
                 ):
        gr.hier_block2.__init__(self,
            "mimo_ofdm_rx_cb",
            gr.io_signature(2, 2, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(1, 1, gr.sizeof_char)) # Output signature

        """
        Parameter initalization
        """
        self.n = n
        self.mimo_technique = mimo_technique
        self.fft_len = fft_len
        self.cp_len = cp_len
        self.frame_length_tag_key = frame_length_tag_key
        self.packet_length_tag_key = packet_length_tag_key
        self.occupied_carriers = occupied_carriers
        self.pilot_carriers =pilot_carriers,
        self.pilot_symbols = pilot_symbols
        self.bps_header = bps_header
        self.bps_payload = bps_payload


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
        Synchronization of the superposition of the N received signals
        """
        add = blocks.add_cc()
        sum_sync_detect = digital.ofdm_sync_sc_cfb(fft_len, cp_len)
        mimo_sync = digital.mimo_ofdm_synchronizer_fbcvc(self.n,
                                                         self.fft_len,
                                                         self.cp_len,
                                                         self.sync_word1,
                                                         self.sync_word2)
        # for i in range(0, self.n): TODO enable
        #     self.connect((self, i), blocks.multiply_const_cc(1.0 / np.sqrt(self.fft_len)), (add, i))
        #     self.connect((self, i), blocks.multiply_const_cc(1.0 / np.sqrt(self.fft_len)), (mimo_sync, 3+i))
        # self.connect(add, sum_sync_detect)
        # self.connect((sum_sync_detect, 0), (mimo_sync, 0))
        # self.connect((sum_sync_detect, 1), (mimo_sync, 1))
        # self.connect(add, (mimo_sync, 2))

        ## Disable sync TODO disable
        dump_cp1 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len + cp_len, cp_len)
        dump_cp2 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len + cp_len, cp_len)
        dump_sync1 = blocks.keep_m_in_n(gr.sizeof_gr_complex * fft_len, 2, 4, 2)
        dump_sync2 = blocks.keep_m_in_n(gr.sizeof_gr_complex * fft_len, 2, 4, 2)
        mult1 = blocks.multiply_const_cc(1.0 / np.sqrt(fft_len))
        mult2 = blocks.multiply_const_cc(1.0 / np.sqrt(fft_len))
        v2s1 = blocks.stream_to_vector(gr.sizeof_gr_complex, self.fft_len)
        v2s2 = blocks.stream_to_vector(gr.sizeof_gr_complex, self.fft_len)
        fft1 = fft.fft_vcc(self.fft_len, True, (), True)
        fft2 = fft.fft_vcc(self.fft_len, True, (), True)
        self.connect((self, 0), dump_cp1, mult1, v2s1, fft1, dump_sync1)
        self.connect((self, 1), dump_cp2, mult2, v2s2, fft2, dump_sync2)


        """
        OFDM demod and MIMO channel estimation
        """
        ofdm_demod = []
        for i in range(0, self.n):
            ofdm_demod.append(fft.fft_vcc(self.fft_len, True, (), True))
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(
            n=self.n,
            fft_len=fft_len,
            pilot_symbols=self.pilot_symbols[:self.n, :self.n],
            pilot_carriers=[-21, -7, 7, 21],
            occupied_carriers=range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22, 27))
        #for i in range(0, self.n): TODO enable
        #    self.connect((mimo_sync, i), ofdm_demod[i], (channel_est, i)) #TODO add coarse freq sync after FFT
        # TODO disable
        stream2tagged_stream_sync1 = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, fft_len, ((48+144)/2)/48, "start")
        stream2tagged_stream_sync2 = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, fft_len, ((48+144)/2)/48, "start")
        self.connect(dump_sync1, stream2tagged_stream_sync1, (channel_est, 0))
        self.connect(dump_sync2, stream2tagged_stream_sync2, (channel_est, 1))
        #self.connect((channel_est, 0), blocks.tag_debug(gr.sizeof_gr_complex*len(_def_occupied_carriers[0]), 'Chanest_rx', 'start'))

        """
        MIMO decoder
        """
        mimo_decoder = mimo_decoder_cc(
            N=self.n,
            mimo_technique=self.mimo_technique,
            vlen=len(range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22, 27)))
        for i in range(0, self.n):
            self.connect((channel_est, i), (mimo_decoder, i))

        """
        Header reader
        """
        header_constellation = _get_constellation(bps_header)
        header_formatter = digital.packet_header_ofdm(
            occupied_carriers, 1,
            packet_length_tag_key,
            frame_length_tag_key,
            packet_num_tag_key,
            bps_header, bps_payload)
        header_reader = digital.mimo_ofdm_header_reader_cc(header_constellation.base(),
                                                           header_formatter.formatter())
        # TODO enable
        self.connect(mimo_decoder, header_reader)
         # TODO disable
        #stream2tagged_stream1 = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 1, 144+48, "start")
        #stream2tagged_stream2 = blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 1, 144, "packet_length")
        #header_ersatz = blocks.keep_m_in_n(gr.sizeof_gr_complex, 144, 144+48, 48)
        #self.connect(mimo_decoder, stream2tagged_stream1, header_reader)
        #self.connect(mimo_decoder, header_ersatz, stream2tagged_stream2)
        #self.connect((channel_est, 0), blocks.tag_debug(gr.sizeof_gr_complex*48, 'Channel est1', 'start'))
        #self.connect((channel_est, 1), blocks.tag_debug(gr.sizeof_gr_complex*48, 'Channel est2', 'start'))

        """
        Payload demod
        """
        payload_constellation = _get_constellation(bps_payload)
        payload_demod = digital.constellation_decoder_cb(payload_constellation.base())
        #self.connect(payload_demod, blocks.tag_debug(gr.sizeof_char, 'demod'))
        payload_pack = blocks.repack_bits_bb(bps_payload, 8, self.packet_length_tag_key, True)
        #self.connect(payload_pack, blocks.tag_debug(gr.sizeof_char, 'pack'))
        crc = digital.crc32_bb(True, self.packet_length_tag_key)
        #self.connect(crc, blocks.tag_debug(gr.sizeof_char, 'crc'))
        self.connect(header_reader, payload_demod, payload_pack, crc, self)
