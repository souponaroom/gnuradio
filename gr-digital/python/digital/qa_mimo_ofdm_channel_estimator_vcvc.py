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

_walsh_sequences = [
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1,-1, 1,-1, 1,-1, 1,-1],
    [1, 1,-1,-1, 1, 1,-1,-1],
    [1,-1,-1, 1, 1,-1,-1, 1],
    [1, 1, 1, 1,-1,-1,-1,-1],
    [1,-1, 1,-1,-1, 1,-1, 1],
    [1, 1,-1,-1,-1,-1, 1, 1],
    [1,-1,-1, 1,-1, 1, 1,-1]
]

class qa_mimo_ofdm_channel_estimator_vcvc (gr_unittest.TestCase):

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
        N=2
        M=2
        channel_matrix = np.array([[1j, -2], [1, 2+2j]])#(np.random.randn(N, M) + 1j * np.random.randn(N, M))
        print 'channel'
        print channel_matrix

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
            m=M, mimo_technique="alamouti"
        )
        static_channel = blocks.multiply_matrix_cc(channel_matrix)
        fft1 = fft.fft_vcc(fft_len, True, (), True)
        fft2 = fft.fft_vcc(fft_len, True, (), True)
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(N, fft_len, _walsh_sequences, [-21, -7, 7, 21])
        dump_cp1 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len+cp_len, cp_len)
        dump_cp2 = blocks.keep_m_in_n(gr.sizeof_gr_complex, fft_len, fft_len + cp_len, cp_len)
        dump_sync1 = blocks.keep_m_in_n(gr.sizeof_gr_complex*fft_len, 6, 8, 2)
        dump_sync2 = blocks.keep_m_in_n(gr.sizeof_gr_complex * fft_len, 6, 8, 2)
        head = blocks.head(gr.sizeof_gr_complex * fft_len, 4)
        sink1 = blocks.vector_sink_c(vlen=fft_len)
        sink2 = blocks.vector_sink_c(vlen=fft_len)

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
        csi = np.empty(shape=[fft_len, N, N], dtype=complex)
        for k in range(0, fft_len):
            for n in range(0, N):
                for m in range(0, M):
                    csi[k][n][m] = pmt.c32vector_ref(pmt.vector_ref(pmt.vector_ref(sink1.tags()[0].value, k), n), m)
        print 'csi'
        print csi


if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_channel_estimator_vcvc, "qa_mimo_ofdm_channel_estimator_vcvc.xml")
