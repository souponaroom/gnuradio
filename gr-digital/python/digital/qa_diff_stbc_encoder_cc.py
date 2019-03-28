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
    def encode(self, basis, input, data_len, block_len, phase_shift):
        mapping_coeffs = np.empty(shape=[2, block_len, data_len/2], dtype=complex)
        output = np.empty(shape=[2, data_len, block_len], dtype=complex)
        # Calculate the coefficients for the new basis.
        for k in range(0, block_len):
            mapping_coeffs[0][k] = input[k::2*block_len] *  np.conj(basis[0]) + input[1*block_len+k::2*block_len] * np.conj(basis[1])
            mapping_coeffs[1][k] = input[k::2*block_len] * -basis[1] + input[1*block_len+k::2*block_len] * basis[0]


        # Calculate the expected result.
        # The first vector is calculated with help of the fixed predecessor.
        pre = M_SQRT_2 * np.exp(1j * phase_shift)
        for k in range(0, block_len):
            output[0][0][k] = mapping_coeffs[0][k][0] * pre - mapping_coeffs[1][k][0] * np.conj(pre)
            output[1][0][k] = mapping_coeffs[0][k][0] * pre + mapping_coeffs[1][k][0] * np.conj(pre)
            output[0][1][k] = -np.conj(output[1][0][k])
            output[1][1][k] = np.conj(output[0][0][k])
        # Iteratively calculate the coefficients for the new vector basis.
        for i in range(1, data_len/2):
            for k in range(0, block_len):
                output[0][i*2][k] = mapping_coeffs[0][k][i]*output[0][(i-1)*2][k] - mapping_coeffs[1][k][i]*np.conj(output[1][(i-1)*2][k])
                output[1][i*2][k] = mapping_coeffs[0][k][i]*output[1][(i-1)*2][k] + mapping_coeffs[1][k][i]*np.conj(output[0][(i-1)*2][k])
                # Calculate the second element of the output sequence after the rules of Alamouti.
                output[0][i*2+1][k] = -np.conj(output[1][i*2][k])
                output[1][i*2+1][k] =  np.conj(output[0][i*2][k])
        return output

    ''' 5 tests validating the correct output of the encoder with random input data, 
        random modulation order and random phase shift. '''
    def test_001_t(self):
        # Define test params.
        data_length = 10
        repetitions = 1

        for i in range(repetitions):
            block_len = 1#np.random.randint(1,9)
            #modulation_order = np.random.randint(1, 4)
            modulation_order = 1
            phase_shift = 0.0#2.0 * np.pi * np.random.randn()
            # Generate random input data.
            data = M_SQRT_2 * np.exp(1j* (2.0*np.random.randint(0, 2**modulation_order, size=data_length*block_len)/(2.0**modulation_order) + phase_shift))
            basis = np.array([M_SQRT_2*np.exp(1j * phase_shift), M_SQRT_2*np.exp(1j * phase_shift)])

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data)
            encoder = digital.diff_stbc_encoder_cc(phase_shift, block_len)
            sink1 = blocks.vector_sink_c()
            sink2 = blocks.vector_sink_c()
            self.tb.connect(src, encoder, sink1)
            self.tb.connect((encoder, 1), sink2)
            # Run flowgraph.
            self.tb.run()
            # Calculate expected result.
            expected_result = self.encode(basis, data, data_length, block_len, phase_shift)
            print 'test'
            print expected_result[0]
            print sink1.data()
            # Check if the expected result equals the actual result.
            self.assertComplexTuplesAlmostEqual(np.reshape(expected_result[0], data_length*block_len), sink1.data(), 4)
            self.assertComplexTuplesAlmostEqual(np.reshape(expected_result[1], data_length*block_len), sink2.data(), 4)

        self.tb.run ()
        # check data

if __name__ == '__main__':
    gr_unittest.run(qa_diff_stbc_encoder_cc, "qa_diff_stbc_encoder_cc.xml")
