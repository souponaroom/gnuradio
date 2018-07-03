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
import pmt

class qa_mimo_channel_estimator_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def produce_training_seq(self, M, length):
        seq = np.empty([M, length], dtype=complex)
        for (m, i), value in np.ndenumerate(seq):
            seq[m][i] = np.exp(-1j*2*np.pi*m*i/length)
        return seq

    def test_001_t (self):
        M = 2
        N = 2
        training_length = 9

        training_sequence = self.produce_training_seq(M, training_length)
        csi = (np.random.randn(N, M) + 1j * np.random.randn(N, M))

        print 'channel'
        print csi
        print 'training sequence'
        print training_sequence

        # Build flowgraph.
        tags = [gr.tag_utils.python_to_tag((0,
                                            pmt.string_to_symbol("pilot"),
                                            pmt.from_long(0),
                                            pmt.from_long(0)))]
        src = [blocks.vector_source_c(data=training_sequence[0],
                                     repeat=False,
                                     tags=tags)]
        sink = [blocks.tag_debug_make(gr.sizeof_gr_complex, "debug1")]
        for m in range(1, M):
            src.append(blocks.vector_source_c(training_sequence[m]))
            sink.append(blocks.null_sink_make(gr.sizeof_gr_complex))
        channel = blocks.multiply_matrix_cc_make(csi)
        estimator = digital.mimo_channel_estimator_cc(M, N, training_sequence)
        for m in range(0, M):
            self.tb.connect(src[m], (channel, m), (estimator, m), sink[m])
        self.tb.run()


        estimation = pmt.to_python(sink[0].current_tags()[0].value)

        # Check if the expected result equals the actual result.
        for n in range(0, N):
            self.assertComplexTuplesAlmostEqual(csi[n], estimation[n], 2)





if __name__ == '__main__':
    gr_unittest.run(qa_mimo_channel_estimator_cc, "qa_mimo_channel_estimator_cc.xml")
