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

class qa_vblast_encoder_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Function which calculates the expected result.
    def encode_vblast(self, input, M):
        output = np.empty(shape=[M, len(input)/M], dtype=complex)
        # Calculate output data for all ports.
        for m in range(0, M):
            output[m] = input[m::M]
        return output

    ''' 
    5 tests validating the correct output of the encoder with random input data
    and random number of antennas M. '''
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5

        # 5 tests validating the correct output of the encoder with random input data.
        for i in range(repetitions):
            # Dice number of transmit antennas.
            num_outputs = np.random.randint(1, 17)
            # Generate random input data.
            data = np.random.randn(data_length*num_outputs) + 1j * np.random.randn(data_length*num_outputs)

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data)
            vblast = digital.vblast_encoder_cc(num_outputs)

            sink = []
            self.tb.connect(src, vblast)
            for m in range(0, num_outputs):
                sink.append(blocks.vector_sink_c())
                self.tb.connect((vblast, m), sink[m])
            # Run flowgraph.
            self.tb.run()

            # Calculate expected result.
            expected_result = self.encode_vblast(data, num_outputs)

            # Check if the expected result equals the actual result.
            for m in range(0, num_outputs):
                self.assertComplexTuplesAlmostEqual(expected_result[m], sink[m].data(), 4)


if __name__ == '__main__':
    gr_unittest.run(qa_vblast_encoder_cc, "qa_vblast_encoder_cc.xml")
