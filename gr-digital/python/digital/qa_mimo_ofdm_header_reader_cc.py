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
from gnuradio.digital.ofdm_txrx import ofdm_tx
import numpy as np
import pmt

class qa_mimo_ofdm_header_reader_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        # Define test params.
        packet_len = 50
        len_tag_key = 'packet_len'
        fft_len = 64
        cp_len = fft_len/4
        pilot_carriers = [-21, -7, 7, 21]
        occupied_carriers = range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22,27)
        N=2
        M=2
        channel_matrix = (np.random.randn(N, M) + 1j * np.random.randn(N, M))
        self.walsh_sequences = np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, -1, 1, -1, 1, -1, 1, -1],
            [1, 1, -1, -1, 1, 1, -1, -1],
            [1, -1, -1, 1, 1, -1, -1, 1],
            [1, 1, 1, 1, -1, -1, -1, -1],
            [1, -1, 1, -1, -1, 1, -1, 1],
            [1, 1, -1, -1, -1, -1, 1, 1],
            [1, -1, -1, 1, -1, 1, 1, -1]
        ])

        src = blocks.vector_source_b(range(packet_len), True, 1, ())
        s2tagged_stream = blocks.stream_to_tagged_stream(gr.sizeof_char, 1,
                                                         packet_len,
                                                         len_tag_key)
        tx = ofdm_tx(
            fft_len=fft_len, cp_len=cp_len,
            packet_length_tag_key=len_tag_key,
            bps_header=1,
            bps_payload=2,
            rolloff=0,
            debug_log=False,
            scramble_bits=False,
            m=M, mimo_technique="vblast"
        )
        static_channel = blocks.multiply_matrix_cc(channel_matrix)
        fft1 = fft.fft_vcc(fft_len, True, (), True)
        fft2 = fft.fft_vcc(fft_len, True, (), True)
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(n=N,
                                                               fft_len=fft_len,
                                                               pilot_symbols=self.walsh_sequences[:N, :N],
                                                               pilot_carriers=pilot_carriers,
                                                               occupied_carriers=occupied_carriers)
        mimo_decoder = digital.vblast_decoder(num_inputs=N, equalizer_type='ZF')

        dump_cp1 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len+cp_len, cp_len)
        dump_cp2 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len + cp_len, cp_len)
        dump_sync1 = blocks.keep_m_in_n(gr.sizeof_gr_complex*fft_len, 6, 8, 2)
        dump_sync2 = blocks.keep_m_in_n(gr.sizeof_gr_complex * fft_len, 6, 8, 2)
        head = blocks.head(gr.sizeof_gr_complex * len(occupied_carriers), 4)
        sink1 = blocks.vector_sink_c(vlen=len(occupied_carriers))
        sink2 = blocks.vector_sink_c(vlen=len(occupied_carriers))

        self.tb.connect(src,
                        s2tagged_stream,
                        tx,
                        (static_channel, 0),
                        dump_cp1,
                        blocks.multiply_const_cc(1.0 / np.sqrt(fft_len)),
                        blocks.stream_to_vector(gr.sizeof_gr_complex, fft_len),
                        fft1,
                        dump_sync1,
                        (channel_est, 0),
                        sink1)
        self.tb.connect((tx, 1),
                        (static_channel, 1),
                        dump_cp2,
                        blocks.multiply_const_cc(1.0 / np.sqrt(fft_len)),
                        blocks.stream_to_vector(gr.sizeof_gr_complex, fft_len),
                        fft2,
                        dump_sync2,
                        (channel_est, 1),
                        head,
                        sink2)
        self.tb.run ()
        # check data



if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_header_reader_cc, "qa_mimo_ofdm_header_reader_cc.xml")
