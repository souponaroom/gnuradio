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

M_SQRT_2 = 1.0/np.sqrt(2)

class qa_diff_stbc_encoder_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Function which encodes the input data to compare it with the actual result of the encoder block.
    def encode(self, basis, input, phase_shift):
        mapping_coeffs = np.empty(shape=[2, len(input)/2], dtype=complex)
        output = np.empty(shape=[2, len(input)], dtype=complex)
        # Calculate the coefficients for the new basis.
        mapping_coeffs[0] = input[::2] *  np.conj(basis[0]) + input[1::2] * np.conj(basis[1])
        mapping_coeffs[1] = input[::2] * -basis[1] + input[1::2] * basis[0]

        # Calculate the expected result.
        # The first vector is calculated with help of the fixed predecessor.
        pre = M_SQRT_2 * np.exp(1j * phase_shift)
        output[0][0] = mapping_coeffs[0][0] * pre - mapping_coeffs[1][0] * np.conj(pre)
        output[1][0] = mapping_coeffs[0][0] * pre + mapping_coeffs[1][0] * np.conj(pre)
        # Iteratively calculate the coefficients for the new vector basis.
        for i in range(1, len(input)/2):
            output[0][i*2] = mapping_coeffs[0][i]*output[0][(i-1)*2] - mapping_coeffs[1][i]*np.conj(output[1][(i-1)*2])
            output[1][i*2] = mapping_coeffs[0][i]*output[1][(i-1)*2] + mapping_coeffs[1][i]*np.conj(output[0][(i-1)*2])
        # Calculate the second element of the output sequence after the rules of Alamouti.
        output[0][1::2] = -np.conj(output[1][::2])
        output[1][1::2] =  np.conj(output[0][::2])
        return output

    ''' 5 tests validating the correct output of the encoder with random input data, 
        random modulation order and random phase shift. '''
    def test_001_t(self):
        # Define test params.
        data_length = 20
        repetitions = 5

        for i in range(repetitions):
            modulation_order = np.random.randint(1, 4)
            phase_shift = 2.0 * np.pi * np.random.randn()
            # Generate random input data.
            data = M_SQRT_2 * np.exp(1j* (2.0*np.pi*np.random.randint(0, 2**modulation_order, size=data_length)/(2.0**modulation_order) + phase_shift))
            basis = np.array([M_SQRT_2*np.exp(1j * phase_shift), M_SQRT_2*np.exp(1j * phase_shift)])

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data)
            encoder = digital.diff_stbc_encoder_cc(phase_shift)
            sink1 = blocks.vector_sink_c()
            sink2 = blocks.vector_sink_c()
            self.tb.connect(src, encoder, sink1)
            self.tb.connect((encoder, 1), sink2)
            # Run flowgraph.
            self.tb.run()

            # Calculate expected result.
            expected_result = self.encode(basis, data, phase_shift)

            # Check if the expected result equals the actual result.
            self.assertComplexTuplesAlmostEqual(expected_result[0], sink1.data(), 4)
            self.assertComplexTuplesAlmostEqual(expected_result[1], sink2.data(), 4)

        self.tb.run ()
        # check data

if __name__ == '__main__':
    gr_unittest.run(qa_diff_stbc_encoder_cc, "qa_diff_stbc_encoder_cc.xml")
