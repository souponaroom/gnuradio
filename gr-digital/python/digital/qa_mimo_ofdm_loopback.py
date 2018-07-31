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

from gnuradio import gr, gr_unittest
from gnuradio import blocks
from gnuradio import fft
import digital_swig as digital
from gnuradio.digital.mimo_encoder_cc import mimo_encoder_cc
from gnuradio.digital.mimo_decoder_cc import mimo_decoder_cc
import numpy as np
import pmt

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
_def_occupied_carriers = (range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22, 27),)
_def_pilot_carriers=((-21, -7, 7, 21,),)
_seq_seed = 42


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
    sw1 = [bpsk[np.random.randint(2)] if x in active_carriers and x % 2 else 0 for x in
           range(fft_len)]
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


class qa_mimo_ofdm_loopback (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def _get_constellation(self, bps):
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


    ''' Generate header and modulate payload data.'''
    def generate_mod_packets(self, data, packet_lengths,
                             occupied_carriers,
                             bps_header, bps_payload,
                             packet_len_tag_key,
                             header_constellation, payload_constellation):
        # Append stream tags to data.
        tags = []
        offset = 0
        for i in range(0, packet_lengths.size):
            tags.append(gr.tag_utils.python_to_tag((offset,
                                                    pmt.string_to_symbol(packet_len_tag_key),
                                                    pmt.from_long(packet_lengths[i]))))
            offset += packet_lengths[i]
        payload_src = blocks.vector_source_b(data=data, tags=tags)

        # Generate header.
        formatter_object_tx = digital.packet_header_ofdm(
            occupied_carriers=occupied_carriers, n_syms=1,
            bits_per_header_sym=bps_header,
            bits_per_payload_sym=bps_payload,
            scramble_header=False)
        header_gen = digital.packet_headergenerator_bb(formatter_object_tx.base(),
                                                       packet_len_tag_key)

        header_mod = digital.chunks_to_symbols_bc(header_constellation.points())
        header_payload_mux = blocks.tagged_stream_mux(
            itemsize=gr.sizeof_gr_complex * 1,
            lengthtagname=packet_len_tag_key,
            tag_preserve_head_pos=1)  # Head tags on the payload stream stay on the head
        self.tb.connect(header_gen, header_mod, (header_payload_mux, 0))

        # Payload modulation
        payload_unpack = blocks.repack_bits_bb(
            8,  # Unpack 8 bits per byte
            bps_payload,
            packet_len_tag_key)

        payload_mod = digital.chunks_to_symbols_bc(payload_constellation.points())

        self.tb.connect(payload_src, payload_unpack, payload_mod, (header_payload_mux, 1))
        self.tb.connect(payload_src, header_gen)
        #self.tb.connect(header_payload_mux, blocks.tag_debug(gr.sizeof_gr_complex, 'Header_Payload_MUX'))
        sink = blocks.vector_sink_c()
        self.tb.connect(header_payload_mux, sink)
        self.tb.run()
        return sink

    def mimo_enc(self, data, tags, mimo_technique, m):
        packet_src = blocks.vector_source_c(data=data, tags=tags)
        mimo_encoder = mimo_encoder_cc(
            M=m,
            mimo_technique=mimo_technique)

        self.tb.connect(packet_src, mimo_encoder)
        #self.tb.connect((mimo_encoder, 0), blocks.tag_debug(gr.sizeof_gr_complex, 'MIMO_encoder 1'))
        #self.tb.connect((mimo_encoder, 1), blocks.tag_debug(gr.sizeof_gr_complex, 'MIMO_encoder 2'))

        sinks = []
        for i in range(0, m):
            sinks.append(blocks.vector_sink_c())
            self.tb.connect((mimo_encoder, i), sinks[i])

        self.tb.run()
        return sinks

    def ofdm_mod(self, data, tags, m, fft_len, cp_len, sync_words, packet_length_tag_key,
                      rolloff=0):
        mimo_sources = []
        for i in range(0, m):
            mimo_sources.append(blocks.vector_source_c(data=data[i], tags=tags[i]))

        allocator = []
        ffter = []
        cyclic_prefixer = []
        normalize = []
        tx_sink = []
        for i in range(0, m):
            mimo_pilot_symbols = np.repeat(_walsh_sequences[i][:m], 4).reshape((m, 4))
            allocator.append(
                digital.ofdm_carrier_allocator_cvc(
                    fft_len,
                    occupied_carriers=_def_occupied_carriers,
                    pilot_carriers=_def_pilot_carriers,
                    pilot_symbols=mimo_pilot_symbols,
                    sync_words=sync_words,
                    len_tag_key=packet_length_tag_key
                )
            )
            ffter.append(
                fft.fft_vcc(
                    fft_len,
                    False,  # Inverse FFT
                    (),  # No window
                    True  # Shift
                )
            )
            cyclic_prefixer.append(
                digital.ofdm_cyclic_prefixer(
                    fft_len,
                    fft_len + cp_len,
                    rolloff,
                    packet_length_tag_key
                )
            )
            normalize.append(blocks.multiply_const_cc(1.0 / np.sqrt(fft_len)))
            self.tb.connect(mimo_sources[i],
                            allocator[i],
                            ffter[i],
                            cyclic_prefixer[i],
                            normalize[i])
            tx_sink.append(blocks.vector_sink_c())
            self.tb.connect(normalize[i], tx_sink[i])
        self.tb.run()
        return tx_sink

    def ofdm_demod(self, data, tags, n, fft_len, cp_len):
        src = []
        dump_cp = []
        dump_sync = []
        mult = []
        v2s = []
        ffter = []
        sinks = []
        for i in range(0, n):
            src.append(blocks.vector_source_c(data=data[i], tags=tags))
            dump_cp.append(blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len + cp_len, cp_len))
            dump_sync.append(blocks.keep_m_in_n(gr.sizeof_gr_complex * fft_len, 1, 3, 2))
            mult.append(blocks.multiply_const_cc(1.0 / np.sqrt(fft_len)))
            v2s.append(blocks.stream_to_vector(gr.sizeof_gr_complex, fft_len))
            ffter.append(fft.fft_vcc(fft_len, True, (), True))
            sinks.append(blocks.vector_sink_c(vlen=fft_len))
            self.tb.connect(src[i], dump_cp[i], mult[i], v2s[i], ffter[i], dump_sync[i], sinks[i])
        self.tb.connect(src[0], blocks.tag_debug(gr.sizeof_gr_complex, 'After sync'))
        self.tb.run()
        return sinks

    def mimo_decode(self, data, tags, n, mimo_technique, fft_len):
        # MIMO channel estimation
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(
            n=n,
            fft_len=fft_len,
            pilot_symbols=_walsh_sequences[:n, :n],
            pilot_carriers=_def_pilot_carriers[0],
            occupied_carriers=_def_occupied_carriers[0])
        # MIMO decoder
        mimo_decoder = mimo_decoder_cc(
            N=n,
            mimo_technique=mimo_technique,
            vlen=len(_def_occupied_carriers[0]))
        src = []
        for i in range(0, n):
            src.append(blocks.vector_source_c(data=data[i], tags=tags[i], vlen=fft_len))
            self.tb.connect(src[i], (channel_est, i), (mimo_decoder, i))
        sink = blocks.vector_sink_c()
        self.tb.connect(mimo_decoder, sink)
        self.tb.run()
        return sink

    '''Accepts tagged stream data and reads/parses the header and demodulates the payload.'''
    def parse(self, data, tags,
              occupied_carriers,
              bps_header, bps_payload,
              header_constellation, payload_constellation
              ):

        rx_src = blocks.vector_source_c(data=data,
                                        tags=tags)
        # Header reader.
        header_formatter_rx = digital.packet_header_ofdm(
            occupied_carriers, 1,
            "packet_length",
            "frame_length",
            "packet_num",
            bps_header, bps_payload)
        header_reader = digital.mimo_ofdm_header_reader_cc(header_constellation.base(),
                                                           header_formatter_rx.formatter())
        # Payload demodulator.
        payload_demod = digital.constellation_decoder_cb(payload_constellation.base())
        payload_pack = blocks.repack_bits_bb(bps_payload, 8, "", True)
        sink_rx = blocks.vector_sink_b()
        self.tb.connect(rx_src, header_reader, payload_demod, payload_pack, sink_rx)
        self.tb.run()
        return sink_rx

    def test_001_t (self):
        # Define test params.
        packet_lengths = np.array([6, 18])
        data_length = np.sum(packet_lengths)
        header_len = len(_def_occupied_carriers[0])
        bps_header = 1
        bps_payload = 1
        fft_len = 64
        cp_len = fft_len/4
        m = 2
        n = 2
        mimo_technique = 'vblast'
        packet_len_tag_key = "packet_length"
        sync_word1 = _make_sync_word1(fft_len, _def_occupied_carriers, _def_pilot_carriers)
        sync_word2 = _make_sync_word2(fft_len, _def_occupied_carriers, _def_pilot_carriers)
        sync_words = [sync_word1, sync_word2]

        # Constellations.
        header_constellation = self._get_constellation(bps_header)
        payload_constellation = self._get_constellation(bps_payload)

        # Payload data.
        data = np.random.randint(0, 100, data_length)

        # Generate modulated packets.
        mod_packets_sink = self.generate_mod_packets(data, packet_lengths,
                                                     _def_occupied_carriers,
                                                     bps_header, bps_payload,
                                                     packet_len_tag_key,
                                                     header_constellation,
                                                     payload_constellation)

        mimo_enc_sinks = self.mimo_enc(data=mod_packets_sink.data(),
                                 tags=mod_packets_sink.tags(),
                                 mimo_technique=mimo_technique,
                                 m=m)

        # pass data
        mimo_enc_data = []
        mimo_enc_tags = []
        for i in range(0, m):
            mimo_enc_data.append(mimo_enc_sinks[i].data())
            mimo_enc_tags.append(mimo_enc_sinks[i].tags())

        mimo_tx = self.ofdm_mod(data=mimo_enc_data,
                                tags=mimo_enc_tags,
                                m=m,
                                fft_len=fft_len,
                                cp_len=cp_len,
                                sync_words=sync_words,
                                packet_length_tag_key=packet_len_tag_key)

        # Channel
        rx_data = []
        for i in range(0, m):
            rx_data.append(mimo_tx[i].data())

        # Simulate tags from sync
        rx_tags = []
        offset = 0
        for i in range(0, packet_lengths.size):
            rx_tags.append(gr.tag_utils.python_to_tag((offset,
                                                       pmt.string_to_symbol("start"),
                                                       pmt.from_long(0))))
            offset += (packet_lengths[i] * 8 / bps_payload + header_len)/n

        # OFDM demod
        ofdm_demod = self.ofdm_demod(data=rx_data,
                                     tags=rx_tags,
                                     n=n,
                                     fft_len=fft_len,
                                     cp_len=cp_len)

        # pass data
        ofdm_demod_data = []
        ofdm_demod_tags = []
        for i in range(0, n):
            ofdm_demod_data.append(ofdm_demod[i].data())
            ofdm_demod_tags.append(ofdm_demod[i].tags())

        # MIMO channel est and decoder
        mimo_dec = self.mimo_decode(data=ofdm_demod_data,
                                    tags=ofdm_demod_tags,
                                    n=n,
                                    mimo_technique=mimo_technique,
                                    fft_len=fft_len)

        # pass data
        mimo_dec_data = mimo_dec.data()
        mimo_dec_tags = mimo_dec.tags()

        # Header and payload reader.
        rx_data = self.parse(mimo_dec_data, mimo_dec_tags,
                             _def_occupied_carriers,
                             bps_header, bps_payload,
                             header_constellation, payload_constellation)


if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_loopback, "qa_mimo_ofdm_loopback.xml")
