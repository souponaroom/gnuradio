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
import pmt
import numpy as np
import random
import cmath

class qa_diversity_combiner_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Test for selection combining with 2 inputs and vector length 1.
    def test_001_t (self):
        mode = 0
        num_inputs = 2
        vlen = 1

        data1 = (np.random.randn(10) + 1j * np.random.randn(10))
        src1 = blocks.vector_source_c(data1)
        data2 = (np.random.randn(10) + 1j * np.random.randn(10))
        src2 = blocks.vector_source_c(data2)
        diversity_combiner = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c()
        self.tb.connect(src1, diversity_combiner, sink)
        self.tb.connect(src2, (diversity_combiner, 1))
        self.tb.run()
        self.assertComplexTuplesAlmostEqual(data1, sink.data(), 6)

    # Test for selection combining with 2 inputs, vector length 1 and random CSI changes.
    def test_003_t (self):
        mode = 0
        num_inputs = 2
        vlen = 1
        num_tags = 4
        tag_pos = [2, 5, 6, 8]

        rand = (np.random.randn(20) + 1j * np.random.randn(20))
        data1 = rand[0:10]
        data2 = rand[10:20]

        tags = []
        expected_result = data1
        for i in range(0, num_tags):
            csi1 = random.random() + 1j*random.random()
            csi2 = random.random() + 1j * random.random()
            csi = pmt.make_c32vector(num_inputs, csi1)
            pmt.c32vector_set(csi, 1, csi2)
            tags.append(gr.tag_utils.python_to_tag((tag_pos[i], pmt.string_to_symbol("csi"), csi, pmt.from_long(0))))

            if(abs(csi1) >= abs(csi2)):
                expected_result[tag_pos[i]:] = data1[tag_pos[i]:]
            else:
                expected_result[tag_pos[i]:] = data2[tag_pos[i]:]

        src1 = blocks.vector_source_c(data=data1,
                                      repeat=False,
                                      vlen=1,
                                      tags=tags)
        src2 = blocks.vector_source_c(data2)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c()
        self.tb.connect(src1, comb, sink)
        self.tb.connect(src2, (comb, 1))
        self.tb.run()
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 6)


if __name__ == '__main__':
    gr_unittest.run(qa_diversity_combiner_cc, "qa_diversity_combiner_cc.xml")
