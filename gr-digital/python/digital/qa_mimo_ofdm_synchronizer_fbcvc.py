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

    def run_flowgraph(self, fft_len, cp_len, symbols, freq_data, trigger_data, ref_data, data1, data2):

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
        head = blocks.head(gr.sizeof_gr_complex * fft_len, symbols)

        self.tb.connect(freq_src, (sync_block, 0))
        self.tb.connect(trigger_src, (sync_block, 1))
        self.tb.connect(ref_src, (sync_block, 2))
        self.tb.connect(src1, (sync_block, 3))
        self.tb.connect(src2, (sync_block, 4))
        self.tb.connect((sync_block, 0), head, sink1)
        self.tb.connect((sync_block, 1), sink2)
        self.tb.run()
        return sink1.data(), sink2.data()

    def rotate_data(self, data, freq, fft_len):
        phase = 0.0
        mod = np.empty(shape=[data.size], dtype=complex)
        for i in range(0, data.size):
            mod[i] = np.exp(1j * phase)
            phase += freq[i]*2.0/fft_len
            phase = phase % (2.0*np.pi)
        return np.multiply(data, mod)

    '''2 tests with no triggers. The data is passed and the cyclic prefix is removed.
    The second test contains a varying fine frequency offset.'''
    def test_001_t (self):
        fft_len = 6
        cp_len = 2
        symbols = 4
        samples = (fft_len+cp_len) * symbols
        ref_data = np.random.randn(samples*2)
        data1 = np.random.randn(samples*2)
        data2 = np.random.randn(samples*2)
        trigger_data = [0]*samples*2
        freq_data = [0.]*samples*2

        # Test without frequency offset.
        res1, res2 = self.run_flowgraph(fft_len=fft_len, cp_len=cp_len,
                                        symbols=symbols, freq_data=freq_data,
                                        trigger_data=trigger_data, ref_data=ref_data,
                                        data1=data1, data2=data2)

        self.assertComplexTuplesAlmostEqual(
            np.delete(data1, np.append(
                np.arange(0, data1.size, fft_len+cp_len),
                np.arange(1, data1.size, fft_len+cp_len))
                      )[:symbols*fft_len], res1[:symbols*fft_len], 2)
        self.assertComplexTuplesAlmostEqual(
            np.delete(data2, np.append(
                np.arange(0, data2.size, fft_len+cp_len),
                np.arange(1, data2.size, fft_len+cp_len))
                      )[:symbols*fft_len], res2[:symbols*fft_len], 2)

        # Test with frequency offset.
        freq_data = np.random.uniform(-3.0, 3.0, samples * 2)
        rotated_data1 = self.rotate_data(data1, freq_data, fft_len)
        rotated_data2 = self.rotate_data(data2, freq_data, fft_len)

        res1, res2 = self.run_flowgraph(fft_len=fft_len, cp_len=cp_len,
                                        symbols=symbols, freq_data=freq_data,
                                        trigger_data=trigger_data, ref_data=ref_data,
                                        data1=rotated_data1, data2=rotated_data2)
        self.assertComplexTuplesAlmostEqual(
            np.delete(data1, np.append(
                np.arange(0, data1.size, fft_len+cp_len),
                np.arange(1, data1.size, fft_len+cp_len))
                      )[:symbols*fft_len], res1[:symbols*fft_len], 2)
        self.assertComplexTuplesAlmostEqual(
            np.delete(data2, np.append(
                np.arange(0, data2.size, fft_len+cp_len),
                np.arange(1, data2.size, fft_len+cp_len))
                      )[:symbols*fft_len], res2[:symbols*fft_len], 2)

    ''' 2 tests with triggers on sync positions (beginning of a symbol). 
    The data is passed and the cyclic prefix is removed.
    After a trigger, the first 2 sync symbols are removed.
    The second test contains a varying fine frequency offset.'''
    def test_002_t (self):
        fft_len = 6
        cp_len = 2
        symbols = 10
        samples = (fft_len+cp_len) * symbols
        ref_data = np.random.randn(samples*2)
        data1 = np.random.randn(samples*2)
        data2 = np.random.randn(samples*2)
        trigger_data = [0]*samples*2
        trigger_data[1*(fft_len+cp_len)] = 1
        trigger_data[5*(fft_len+cp_len)] = 1
        freq_data = [0.]*samples*2

        # Test without freqency offset.
        res1, res2 = self.run_flowgraph(fft_len=fft_len, cp_len=cp_len,
                                        symbols=symbols, freq_data=freq_data,
                                        trigger_data=trigger_data, ref_data=ref_data,
                                        data1=data1, data2=data2)

        expected_result1 = np.delete(np.delete(data1, np.append(
            np.arange(0, data1.size, fft_len+cp_len),
            np.arange(1, data1.size, fft_len+cp_len))),
                                     np.append(np.arange(fft_len, 3*fft_len), np.arange(5*fft_len, 7*fft_len)))
        expected_result2 = np.delete(np.delete(data2, np.append(
            np.arange(0, data2.size, fft_len + cp_len),
            np.arange(1, data2.size, fft_len + cp_len))),
                                     np.append(np.arange(fft_len, 3*fft_len), np.arange(5*fft_len, 7*fft_len)))
        self.assertComplexTuplesAlmostEqual(
            expected_result1[:fft_len*symbols], res1[:fft_len*symbols], 2)
        self.assertComplexTuplesAlmostEqual(
            expected_result2[:fft_len*symbols], res2[:fft_len*symbols], 2)

        # Test with fine frequency offset.
        freq_data = np.random.uniform(-3.0, 3.0, samples * 2)
        rotated_data1 = self.rotate_data(data1, freq_data, fft_len)
        rotated_data2 = self.rotate_data(data2, freq_data, fft_len)

        res1, res2 = self.run_flowgraph(fft_len=fft_len, cp_len=cp_len,
                                        symbols=symbols, freq_data=freq_data,
                                        trigger_data=trigger_data, ref_data=ref_data,
                                        data1=rotated_data1, data2=rotated_data2)

        expected_result1 = np.delete(np.delete(data1, np.append(
            np.arange(0, data1.size, fft_len+cp_len),
            np.arange(1, data1.size, fft_len+cp_len))),
                                     np.append(np.arange(fft_len, 3*fft_len), np.arange(5*fft_len, 7*fft_len)))
        expected_result2 = np.delete(np.delete(data2, np.append(
            np.arange(0, data2.size, fft_len + cp_len),
            np.arange(1, data2.size, fft_len + cp_len))),
                                     np.append(np.arange(fft_len, 3*fft_len), np.arange(5*fft_len, 7*fft_len)))
        self.assertComplexTuplesAlmostEqual(
            expected_result1[:fft_len*symbols], res1[:fft_len*symbols], 2)
        self.assertComplexTuplesAlmostEqual(
            expected_result2[:fft_len*symbols], res2[:fft_len*symbols], 2)

    '''2 tests with triggers on async positions (not at the beginning of a symbol). 
    The data is passed and the cyclic prefix is removed.
    After a trigger in the middle of a symbol, the symbol is filled from there on with zeros.
    The first 2 sync symbols (in sync of the new trigger) are removed.
    The second test contains a varying fine frequency offset.'''
    def test_003_t (self):
        offset = 3
        fft_len = 6
        cp_len = 2
        symbols = 4
        samples = (fft_len+cp_len) * symbols
        ref_data = np.random.randn(samples*2)
        data1 = np.random.randn(samples*2)
        data2 = np.random.randn(samples*2)
        trigger_data = [0]*samples*2
        trigger_data[offset] = 1
        freq_data = [0.]*samples*2

        # Test without freq offset.
        res1, res2 = self.run_flowgraph(fft_len=fft_len, cp_len=cp_len,
                                        symbols=symbols, freq_data=freq_data,
                                        trigger_data=trigger_data, ref_data=ref_data,
                                        data1=data1, data2=data2)
        exp1_with_cp = np.append(data1[:offset],
                                 np.append(np.zeros(fft_len-offset+cp_len),
                                           data1[offset+2*(fft_len+cp_len):(symbols+2)*(fft_len*cp_len)]))
        exp1 = np.delete(exp1_with_cp,
                         np.append(np.arange(0, data1.size, fft_len+cp_len),
                                   np.arange(1, data1.size, fft_len+cp_len)))
        exp2_with_cp = np.append(data2[:offset],
                                 np.append(np.zeros(fft_len - offset + cp_len),
                                           data2[offset + 2 * (fft_len + cp_len):(symbols + 2)*(fft_len*cp_len)]))
        exp2 = np.delete(exp2_with_cp,
                         np.append(np.arange(0, data2.size, fft_len + cp_len),
                                   np.arange(1, data2.size, fft_len + cp_len)))

        self.assertComplexTuplesAlmostEqual(exp1[:fft_len*symbols], res1[:fft_len*symbols], 2)
        self.assertComplexTuplesAlmostEqual(exp2[:fft_len*symbols], res2[:fft_len * symbols], 2)

        # Test with fine frequency offset.
        freq_data = np.random.uniform(-3.0, 3.0, samples * 2)
        rotated_data1 = self.rotate_data(data1, freq_data, fft_len)
        rotated_data2 = self.rotate_data(data2, freq_data, fft_len)

        res1, res2 = self.run_flowgraph(fft_len=fft_len, cp_len=cp_len,
                                        symbols=symbols, freq_data=freq_data,
                                        trigger_data=trigger_data, ref_data=ref_data,
                                        data1=rotated_data1, data2=rotated_data2)
        exp1_with_cp = np.append(data1[:offset],
                                 np.append(np.zeros(fft_len-offset+cp_len),
                                           data1[offset+2*(fft_len+cp_len):(symbols+2)*(fft_len*cp_len)]))
        exp1 = np.delete(exp1_with_cp,
                         np.append(np.arange(0, data1.size, fft_len+cp_len),
                                   np.arange(1, data1.size, fft_len+cp_len)))
        exp2_with_cp = np.append(data2[:offset],
                                 np.append(np.zeros(fft_len - offset + cp_len),
                                           data2[offset + 2 * (fft_len + cp_len):(symbols + 2)*(fft_len*cp_len)]))
        exp2 = np.delete(exp2_with_cp,
                         np.append(np.arange(0, data2.size, fft_len + cp_len),
                                   np.arange(1, data2.size, fft_len + cp_len)))

        self.assertComplexTuplesAlmostEqual(exp1[:fft_len*symbols], res1[:fft_len*symbols], 2)
        self.assertComplexTuplesAlmostEqual(exp2[:fft_len*symbols], res2[:fft_len * symbols], 2)


if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_synchronizer_fbcvc, "qa_mimo_ofdm_synchronizer_fbcvc.xml")
