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

class qa_diversity_combiner_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Function that dices the csi vectors and calculates the expected result.
    def dice_csi_tags (self, num_inputs, data, num_tags, tag_pos):
        tags = []
        expected_result = data[0]
        for i in range(0, num_tags):
            # Dice CSI for one symbol.
            csi = (np.random.randn(num_inputs) + 1j * np.random.randn(num_inputs))
            # Assign the CSI vector to a PMT vector.
            csi_pmt = pmt.make_c32vector(num_inputs, csi[0])
            for j, channel in enumerate(csi):
                pmt.c32vector_set(csi_pmt, j, channel)
            # Append stream tags with CSI to data stream.
            tags.append(gr.tag_utils.python_to_tag((tag_pos[i],
                                                    pmt.string_to_symbol("csi"),
                                                    csi_pmt,
                                                    pmt.from_long(0))))
            # Calculate the expected result by selecting the channel with the greater magnitude.
            expected_result[tag_pos[i]:] = data[np.argmax(np.abs(csi))][tag_pos[i]:]
        return tags, expected_result

    # Test for selection combining with 2 inputs, vector length 1 and 4 random CSI changes.
    def test_001_t (self):
        # Define test params.
        mode = 0
        num_inputs = 2
        vlen = 1
        tag_pos = [2, 5, 6, 8]
        num_tags = len(tag_pos)

        data1 = (np.random.randn(10*vlen) + 1j * np.random.randn(10*vlen))
        data2 = (np.random.randn(10*vlen) + 1j * np.random.randn(10*vlen))
        data = [data1, data2]

        tags, expected_result = self.dice_csi_tags(num_inputs, data, num_tags, tag_pos)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data1,
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        src2 = blocks.vector_source_c(data=data2)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c()
        self.tb.connect(src1, comb, sink)
        self.tb.connect(src2, (comb, 1))
        self.tb.run()
        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 6)

    # Test for selection combining with 3 inputs, vector length 1 and 5 random CSI changes.
    def test_002_t (self):
        # Define test params.
        mode = 0
        num_inputs = 3
        vlen = 1
        tag_pos = [1, 3, 6, 7, 8]
        num_tags = len(tag_pos)

        data1 = (np.random.randn(10*vlen) + 1j * np.random.randn(10*vlen))
        data2 = (np.random.randn(10*vlen) + 1j * np.random.randn(10*vlen))
        data3 = (np.random.randn(10*vlen) + 1j * np.random.randn(10*vlen))
        data = [data1, data2, data3]

        tags, expected_result = self.dice_csi_tags(num_inputs, data, num_tags, tag_pos)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data1,
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        src2 = blocks.vector_source_c(data=data2,
                                      vlen=vlen)
        src3 = blocks.vector_source_c(data=data3,
                                      vlen=vlen)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c(vlen=vlen)
        self.tb.connect(src1, comb, sink)
        self.tb.connect(src2, (comb, 1))
        self.tb.connect(src3, (comb, 2))
        self.tb.run()
        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 6)


if __name__ == '__main__':
    gr_unittest.run(qa_diversity_combiner_cc, "qa_diversity_combiner_cc.xml")
