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

    # Function which randomly generates the CSI vectors and calculates the expected result.
    def dice_csi_tags (self, num_inputs, data, vlen, num_tags, tag_pos, combining_technique):
        tags = []
        if combining_technique == 'SC':
            # SC initially selects channel 0 before it is updated by the CSI.
            expected_result = np.array(data[0], copy=True)
        else:
            # MRC initially weights all channels equally before it is updated by the CSI.
            expected_result = np.dot(np.full(num_inputs, 1.0/num_inputs), data)

        for i in range(0, num_tags):
            # Randomly generate CSI for one symbol.
            csi = (np.random.randn(num_inputs) + 1j * np.random.randn(num_inputs))
            # Assign the CSI vector to a PMT vector.
            csi_pmt = pmt.make_c32vector(num_inputs, csi[0])
            for k, channel in enumerate(csi):
                pmt.c32vector_set(v=csi_pmt, k=k, x=channel)
            # Append stream tags with CSI to data stream.
            tags.append(gr.tag_utils.python_to_tag((tag_pos[i],
                                                    pmt.string_to_symbol("csi"),
                                                    csi_pmt,
                                                    pmt.from_long(0))))
            # Calculate the expected result.
            if combining_technique == 'SC':
                # Select the path with the greatest magnitude.
                expected_result[vlen*tag_pos[i]:] = data[np.argmax(np.abs(csi))][vlen*tag_pos[i]:]
            elif combining_technique == 'MRC':
                # Calculate the normalized weighting vector.
                weighting_vector = np.sqrt(np.divide(np.square(np.abs(csi)),
                                                     np.sum(np.square(np.abs(csi))))) * np.exp(-1j * np.angle(csi))
                # The weighted sum corresponds to a matrix multiplication
                # of the weighting vector and the data matrix.
                expected_result[vlen*tag_pos[i]:] = (np.dot(weighting_vector, data))[vlen*tag_pos[i]:]
        return tags, expected_result

    # Test for SC: inputs: 2, vector length: 1, random CSI changes: 4.
    def test_001_t (self):
        # Define test params.
        data_length = 20
        mode = 'SC'
        num_inputs = 2
        vlen = 1
        tag_pos = np.array([2, 5, 6, 8])  # Vector-wise indexing.
        num_tags = len(tag_pos)
        # Generate random input data.
        data = np.array(
            [(np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen))])
        for i in range(1, num_inputs):
            data = np.vstack((data, [np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen)]))

        # Generate the CSI vectors and calculate the expected result.
        tags, expected_result = self.dice_csi_tags(num_inputs=num_inputs,
                                                   data=data,
                                                   vlen=vlen,
                                                   num_tags=num_tags,
                                                   tag_pos=tag_pos,
                                                   combining_technique=mode)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data[0],
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c()
        self.tb.connect(src1, comb, sink)
        # Connect all other sources.
        for i in range(1, num_inputs):
            self.tb.connect(blocks.vector_source_c(data=data[i], vlen=vlen), (comb, i))
        # Run flowgraph.
        self.tb.run()
        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

    # Test for SC: inputs: 3, vector length: 2, random CSI changes: 5.
    def test_002_t (self):
        # Define test params.
        data_length = 20
        mode = 'SC'
        num_inputs = 3
        vlen = 2
        tag_pos = [0, 3, 6, 7, 8]  # Vector-wise indexing.
        num_tags = len(tag_pos)
        # Generate random input data.
        data = np.array(
            [(np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen))])
        for i in range(1, num_inputs):
            data = np.vstack((data, [np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen)]))

        # Generate the CSI vectors and calculate the expected result.
        tags, expected_result = self.dice_csi_tags(num_inputs=num_inputs,
                                                   data=data,
                                                   vlen=vlen,
                                                   num_tags=num_tags,
                                                   tag_pos=tag_pos,
                                                   combining_technique=mode)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data[0],
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c(vlen=vlen)
        self.tb.connect(src1, comb, sink)
        # Connect all other sources.
        for i in range(1, num_inputs):
            self.tb.connect(blocks.vector_source_c(data=data[i], vlen=vlen), (comb, i))
        # Run flowgraph.
        self.tb.run()
        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

    # Test for SC: inputs: 8, vector length: 4, random CSI changes: 4.
    def test_003_t(self):
        # Define test params.
        data_length = 20
        mode = 'SC'
        num_inputs = 8
        vlen = 4
        tag_pos = np.array([3, 5, 6, 8])  # Vector-wise indexing.
        num_tags = len(tag_pos)
        # Generate random input data.
        data = np.array([(np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen))])
        for i in range(1, num_inputs):
            data = np.vstack((data, [np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen)]))

        # Generate the CSI vectors and calculate the expected result.
        tags, expected_result = self.dice_csi_tags(num_inputs=num_inputs,
                                                   data=data,
                                                   vlen=vlen,
                                                   num_tags=num_tags,
                                                   tag_pos=tag_pos,
                                                   combining_technique=mode)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data[0],
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c(vlen=vlen)
        self.tb.connect(src1, comb, sink)
        # Connect all other
        for i in range(1, num_inputs):
            self.tb.connect(blocks.vector_source_c(data=data[i], vlen=vlen), (comb, i))
        # Run flowgraph.
        self.tb.run()
        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

    # Test for MRC: inputs: 2, vector length: 1, 4 random CSI changes.
    def test_004_t(self):
        # Define test params.
        data_length = 20
        mode = 'MRC'
        num_inputs = 2
        vlen = 1
        tag_pos = [2, 4, 5, 9]  # Vector-wise indexing.
        num_tags = len(tag_pos)
        # Generate random input data.
        data = np.array(
            [(np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen))])
        for i in range(1, num_inputs):
            data = np.vstack((data, [np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen)]))

        # Generate the CSI vectors and calculate the expected result.
        tags, expected_result = self.dice_csi_tags(num_inputs=num_inputs,
                                                   data=data,
                                                   vlen=vlen,
                                                   num_tags=num_tags,
                                                   tag_pos=tag_pos,
                                                   combining_technique=mode)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data[0],
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c(vlen=vlen)
        self.tb.connect(src1, comb, sink)
        # Connect all other sources.
        for i in range(1, num_inputs):
            self.tb.connect(blocks.vector_source_c(data=data[i], vlen=vlen), (comb, i))
        # Run flowgraph.
        self.tb.run()

        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

    # Test for MRC: inputs: 3, vector length: 2, 6 random CSI changes.
    def test_005_t(self):
        # Define test params.
        data_length = 20
        mode = 'MRC'
        num_inputs = 3
        vlen = 2
        tag_pos = [1, 3, 4, 8, 9]  # Vector-wise indexing.
        num_tags = len(tag_pos)
        # Generate random input data.
        data = np.array(
            [(np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen))])
        for i in range(1, num_inputs):
            data = np.vstack((data, [np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen)]))

        # Generate the CSI vectors and calculate the expected result.
        tags, expected_result = self.dice_csi_tags(num_inputs=num_inputs,
                                                   data=data,
                                                   vlen=vlen,
                                                   num_tags=num_tags,
                                                   tag_pos=tag_pos,
                                                   combining_technique=mode)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data[0],
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c(vlen=vlen)
        self.tb.connect(src1, comb, sink)
        # Connect all other sources.
        for i in range(1, num_inputs):
            self.tb.connect(blocks.vector_source_c(data=data[i], vlen=vlen), (comb, i))
        # Run flowgraph.
        self.tb.run()

        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

    # Test for MRC: inputs: 8, vector length: 5, 2 random CSI changes.
    def test_006_t(self):
        # Define test params.
        data_length = 20
        mode = 'MRC'
        num_inputs = 8
        vlen = 5
        tag_pos = [1, 2]  # Vector-wise indexing.
        num_tags = len(tag_pos)
        # Generate random input data.
        data = np.array(
            [(np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen))])
        for i in range(1, num_inputs):
            data = np.vstack((data, [np.random.randn(data_length * vlen) + 1j * np.random.randn(data_length * vlen)]))

        # Generate the CSI vectors and calculate the expected result.
        tags, expected_result = self.dice_csi_tags(num_inputs=num_inputs,
                                                   data=data,
                                                   vlen=vlen,
                                                   num_tags=num_tags,
                                                   tag_pos=tag_pos,
                                                   combining_technique=mode)

        # Build up the test flowgraph and run it.
        src1 = blocks.vector_source_c(data=data[0],
                                      repeat=False,
                                      vlen=vlen,
                                      tags=tags)
        comb = digital.diversity_combiner_cc(num_inputs, vlen, mode)
        sink = blocks.vector_sink_c(vlen=vlen)
        self.tb.connect(src1, comb, sink)
        # Connect all other sources.
        for i in range(1, num_inputs):
            self.tb.connect(blocks.vector_source_c(data=data[i], vlen=vlen), (comb, i))
        # Run flowgraph.
        self.tb.run()

        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

if __name__ == '__main__':
    gr_unittest.run(qa_diversity_combiner_cc, "qa_diversity_combiner_cc.xml")
