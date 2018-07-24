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
import digital_swig as digital
import numpy as np

class qa_mimo_ofdm_synchronizer_fbcvc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def run_flowgraph(self, fft_len, cp_len, freq_data, trigger_data, ref_data, data1, data2):

        sync_symbol1 = np.random.randn(fft_len)
        sync_symbol2 = np.random.randn(fft_len)

        ref_src = blocks.vector_source_c(ref_data)
        src1 = blocks.vector_source_c(data1)
        src2 = blocks.vector_source_c(data2)
        trigger_src = blocks.vector_source_b(trigger_data)
        freq_src = blocks.vector_source_f(freq_data)

        sync_block = digital.mimo_ofdm_synchronizer_fbcvc(n=2,
                                                          fft_len=fft_len,
                                                          cp_len=cp_len,
                                                          sync_symbol1=sync_symbol1,
                                                          sync_symbol2=sync_symbol2)
        sink1 = blocks.vector_sink_c(vlen=fft_len)
        sink2 = blocks.vector_sink_c(vlen=fft_len)
        head = blocks.head(gr.sizeof_gr_complex * fft_len, 4)

        self.tb.connect(freq_src, (sync_block, 0))
        self.tb.connect(trigger_src, (sync_block, 1))
        self.tb.connect(ref_src, (sync_block, 2))
        self.tb.connect(src1, (sync_block, 3))
        self.tb.connect(src2, (sync_block, 4))
        self.tb.connect((sync_block, 0), head, sink1)
        self.tb.connect((sync_block, 1), sink2)
        self.tb.run()
        return sink1.data(), sink2.data()


    def test_001_t (self):
        fft_len = 6
        cp_len = 2
        length = fft_len * 10
        ref_data = np.random.randn(length)
        data1 = np.random.randn(length)
        data2 = np.random.randn(length)
        trigger_data = [0]*length
        freq_data = [0.]*length

        res1, res2 = self.run_flowgraph(fft_len=fft_len,
                                        cp_len=cp_len,
                                        freq_data=freq_data,
                                        trigger_data=trigger_data,
                                        ref_data=ref_data,
                                        data1=data1,
                                        data2=data2)

        print 'data1'
        print data1
        print 'res1'
        print res1




if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_synchronizer_fbcvc, "qa_mimo_ofdm_synchronizer_fbcvc.xml")
