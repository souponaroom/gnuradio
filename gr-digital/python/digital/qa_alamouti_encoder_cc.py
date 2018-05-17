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

class qa_alamouti_encoder_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Function which calculates the expected result.
    def encode_alamouti(self, input):
        output = np.empty(shape=[2, len(input)], dtype=complex)
        # Calculate output data for port 0.
        output[0][::2] = input[::2]
        output[0][1::2] = -np.conj(input[1::2])
        # Calculate output data for port 1.
        output[1][::2] = input[1::2]
        output[1][1::2] = np.conj(input[::2])

        return output

    # 5 tests validating the correct output of the encoder with random input data.
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5

        for i in range(repetitions):
            # Generate random input data.
            data = np.random.randn(data_length) + 1j * np.random.randn(data_length)

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data)
            alamouti = digital.alamouti_encoder_cc()
            sink1 = blocks.vector_sink_c()
            sink2 = blocks.vector_sink_c()
            self.tb.connect(src, alamouti, sink1)
            self.tb.connect((alamouti, 1), sink2)
            # Run flowgraph.
            self.tb.run()

            # Calculate expected result.
            expected_result = self.encode_alamouti(data)

            # Check if the expected result equals the actual result.
            self.assertComplexTuplesAlmostEqual(expected_result[0], sink1.data(), 4)
            self.assertComplexTuplesAlmostEqual(expected_result[1], sink2.data(), 4)

if __name__ == '__main__':
    gr_unittest.run(qa_alamouti_encoder_cc, "qa_alamouti_encoder_cc.xml")
