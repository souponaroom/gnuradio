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
from gnuradio import analog
from gnuradio import fft
import digital_swig as digital
from gnuradio.digital.ofdm_txrx import ofdm_tx
from gnuradio.digital.mimo_ofdm_rx_cb import mimo_ofdm_rx_cb
import numpy as np

class qa_mimo_ofdm_rx_cb (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None


    def test_001_t (self):
        # Define test params.
        packet_len = 8
        len_tag_key = 'packet_length'
        fft_len = 64
        cp_len = fft_len/4
        N=2
        M=2
        channel_matrix = (np.random.randn(N, M) + 1j * np.random.randn(N, M))
        for i in range(0, 1):

            src = blocks.vector_source_b(range(packet_len*4), True, 1, ())
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
            const = analog.sig_source_c(fft_len, analog.GR_COS_WAVE, 5.0, 1.0)
            mult1 = blocks.multiply_cc()
            mult2 = blocks.multiply_cc()

            rx = mimo_ofdm_rx_cb(
                n=N,
                mimo_technique='vblast',
                fft_len=fft_len,
                cp_len=cp_len,
                packet_length_tag_key=len_tag_key,
                bps_header=1,
                bps_payload=1
            )

            sink = blocks.vector_sink_b()

            self.tb.connect(src, blocks.head(gr.sizeof_char, packet_len*100), s2tagged_stream, tx)
            self.tb.connect((tx, 0), (static_channel, 0), mult1, (rx, 0))
            self.tb.connect((tx, 1), (static_channel, 1), mult2, (rx, 1))
            self.tb.connect(const, (mult1, 1))
            self.tb.connect(const, (mult2, 1))
            # self.tb.connect((tx, 0), chan_sink1)
            # self.tb.connect((tx, 1), chan_sink2)
            # self.tb.connect(chan_src1, (rx, 0))
            # self.tb.connect(chan_src2, (rx, 1))
            self.tb.connect(rx, blocks.head(gr.sizeof_char, (packet_len+4)*4), sink)

            self.tb.run ()
            # check data
            print 'result'
            for i in range(0, len(sink.data()) / (packet_len+4)):
                print sink.data()[i * (packet_len+4):(i + 1) * (packet_len+4)]
            #self.assertComplexTuplesAlmostEqual(range(packet_len*10), sink.data(), 2)



if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_rx_cb, "qa_mimo_ofdm_rx_cb.xml")
