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

class qa_mimo_ofdm_channel_estimator_vcvc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        # Define test params.
        packet_len = 14
        len_tag_key = 'packet_length'
        fft_len = 64
        cp_len = fft_len/4
        pilot_carriers = [-21, -7, 7, 21]
        occupied_carriers = range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22,27)
        N=2
        M=2
        channel_matrix = np.array([[1,0],[0,1]])#(np.random.randn(N, M) + 1j * np.random.randn(N, M))
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

        src = blocks.vector_source_b(range(packet_len*10), False, 1, ())
        s2tagged_stream = blocks.stream_to_tagged_stream(gr.sizeof_char, 1,
                                                         packet_len,
                                                         len_tag_key)
        tx = ofdm_tx(
            fft_len=fft_len, cp_len=cp_len,
            packet_length_tag_key=len_tag_key,
            bps_header=1,
            bps_payload=1,
            rolloff=0,
            debug_log=False,
            scramble_bits=False,
            m=M, mimo_technique="vblast"
        )
        static_channel = blocks.multiply_matrix_cc(channel_matrix)
        print 'channel matrix'
        print channel_matrix
        fft1 = fft.fft_vcc(fft_len, True, (), True)
        fft2 = fft.fft_vcc(fft_len, True, (), True)
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(n=N,
                                                               fft_len=fft_len,
                                                               pilot_symbols=self.walsh_sequences[:N, :N],
                                                               pilot_carriers=pilot_carriers,
                                                               occupied_carriers=occupied_carriers)
        dump_cp1 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len+cp_len, cp_len)
        dump_cp2 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len + cp_len, cp_len)
        dump_sync1 = blocks.keep_m_in_n(gr.sizeof_gr_complex*fft_len, 2, 4, 2)
        dump_sync2 = blocks.keep_m_in_n(gr.sizeof_gr_complex * fft_len, 2, 4, 2)
        #head = blocks.head(gr.sizeof_gr_complex * len(occupied_carriers), 10)
        sink1 = blocks.vector_sink_c(vlen=len(occupied_carriers))
        sink2 = blocks.vector_sink_c(vlen=len(occupied_carriers))
        debug_sink1 = blocks.vector_sink_c(vlen=fft_len)
        debug_sink2 = blocks.vector_sink_c(vlen=fft_len)

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
        #self.tb.connect((channel_est, 0), blocks.tag_debug(gr.sizeof_gr_complex*len(occupied_carriers), 'channel est pls'),)
        self.tb.connect((tx, 1),
                        (static_channel, 1),
                        dump_cp2,
                        blocks.multiply_const_cc(1.0 / np.sqrt(fft_len)),
                        blocks.stream_to_vector(gr.sizeof_gr_complex, fft_len),
                        fft2,
                        dump_sync2,
                        (channel_est, 1),
                        #head,
                        sink2)
        self.tb.connect(dump_sync1, debug_sink1)
        self.tb.connect(dump_sync2, debug_sink2)
        self.tb.run ()
        # check data
        csi = np.empty(shape=[len(occupied_carriers), N, N], dtype=complex)
        for k in range(0, len(occupied_carriers)):
            for n in range(0, N):
                for m in range(0, M):
                    csi[k][n][m] = pmt.c32vector_ref(pmt.vector_ref(pmt.vector_ref(sink1.tags()[0].value, k), n), m)
        print 'csi'
        print csi


        for c in range(0, len(occupied_carriers)):
            for n in range(0, N):
                self.assertComplexTuplesAlmostEqual(channel_matrix[n], csi[c][n], 1)


    # def test_002_t(self):
    #     N=2
    #     fft_len = 4
    #     length = 10
    #     occupied_carriers = [-1]
    #     data = np.random.randn(N, fft_len*10)
    #     src1 = blocks.vector_source_c(data[0], vlen=fft_len)
    #     src2 = blocks.vector_source_c(data[1], vlen=fft_len)
    #     channel_est = digital.mimo_ofdm_channel_estimator_vcvc(n=N,
    #                                                            fft_len=fft_len,
    #                                                            pilot_symbols=[[1,1],[1,-1]],
    #                                                            pilot_carriers=[0],
    #                                                            occupied_carriers=occupied_carriers)
    #     sink1 = blocks.vector_sink_c(vlen=len(occupied_carriers))
    #     sink2 = blocks.vector_sink_c(vlen=len(occupied_carriers))
    #     head = blocks.head(gr.sizeof_gr_complex * len(occupied_carriers), 5)
    #     self.tb.connect(src1, (channel_est, 0), head, sink1)
    #     self.tb.connect(src2, (channel_est, 1), sink2)
    #     self.tb.run()
    #     self.assertComplexTuplesAlmostEqual(data[0][-1 + fft_len / 2:5 * fft_len:fft_len],
    #         sink1.data()[0:len(occupied_carriers) * length / 2], 2)
    #     self.assertComplexTuplesAlmostEqual(data[1][-1 + fft_len / 2:5 * fft_len:fft_len],
    #         sink2.data()[0:len(occupied_carriers) * length / 2], 2)

if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_channel_estimator_vcvc, "qa_mimo_ofdm_channel_estimator_vcvc.xml")
