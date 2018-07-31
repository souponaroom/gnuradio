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

class qa_alamouti_decoder_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Function which randomly generates the CSI vectors and calculates the expected result.
    def dice_csi_tags(self, input, num_tags, tag_pos, vlen=1):
        tags = []
        # Calculate initial behaviour before first tag.
        output = np.empty(shape=[len(input), len(input[0])], dtype=complex)
        for k in range(0, vlen):
            output[k][::2] = (input[k][::2] + np.conj(input[k][1::2]))/2.0
            output[k][1::2] = (input[k][::2] - np.conj(input[k][1::2]))/2.0

        # Iterate over tags and update the calculated output according to the diced CSI.
        for i in range(0, num_tags):
            # Randomly generate CSI for one symbol.
            csi = (np.random.randn(vlen, 1, 2) + 1j * np.random.randn(vlen, 1, 2))
            # Assign the CSI vector to a PMT vector.
            csi_pmt = pmt.make_vector(vlen, pmt.make_vector(1, pmt.make_c32vector(2, 1.0)))
            for k, carrier in enumerate(csi):
                carrier_vector_pmt = pmt.make_vector(1, pmt.make_c32vector(2, csi[k][0][0]))
                for l, rx in enumerate(csi[k]):
                    row_vector_pmt = pmt.make_c32vector(2, csi[k][l][0])
                    for m, tx in enumerate(csi[k][l]):
                        pmt.c32vector_set(v=row_vector_pmt, k=m, x=csi[k][l][m])
                    pmt.vector_set(carrier_vector_pmt, l, row_vector_pmt)
                pmt.vector_set(csi_pmt, k, carrier_vector_pmt)
            # Append stream tags with CSI to data stream.
            tags.append(gr.tag_utils.python_to_tag((tag_pos[i],
                                                    pmt.string_to_symbol("csi"),
                                                    csi_pmt,
                                                    pmt.from_long(0))))
            # Calculate the expected result.
            for k in range(0, vlen):
                total_branch_energy = np.sum(np.square(np.abs(csi[k][0])))
                if tag_pos[i]%2 == 0:
                    output[k][tag_pos[i]  ::2] = (np.conj(csi[k][0][0])*input[k][tag_pos[i]::2] +
                                               csi[k][0][1]*np.conj(input[k][tag_pos[i]+1::2]))/total_branch_energy
                    output[k][tag_pos[i]+1::2] = (np.conj(csi[k][0][1])*input[k][tag_pos[i]::2] -
                                               csi[k][0][0]*np.conj(input[k][tag_pos[i]+1::2]))/total_branch_energy
                else:
                    output[k][tag_pos[i]+1::2] = (np.conj(csi[k][0][0])*input[k][tag_pos[i]+1::2] +
                                               csi[k][0][1]*np.conj(input[k][tag_pos[i]+2::2]))/total_branch_energy
                    output[k][tag_pos[i]+2::2] = (np.conj(csi[k][0][1])*input[k][tag_pos[i]+1::2] -
                                               csi[k][0][0]*np.conj(input[k][tag_pos[i]+2::2]))/total_branch_energy

        return tags, np.reshape(output.T, vlen*len(input[0]))

    ''' 
    5 tests validating the correct output of the decoder with random input data and vector length 1.
    '''
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5
        num_tags = 4
        vlen = 1

        for i in range(repetitions):
            # Generate random input data.
            data = np.random.randn(vlen, data_length) + 1j * np.random.randn(vlen, data_length)
            # Generate random tag positions.
            tag_pos = np.random.randint(low=0, high=data_length/2, size=num_tags)*2
            tag_pos = np.sort(tag_pos)
            # Calculate expected result.
            tags, expected_result = self.dice_csi_tags(data,
                                                       num_tags,
                                                       tag_pos)

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=np.reshape(data.T, vlen*data_length),
                                         repeat=False,
                                         tags=tags)
            alamouti = digital.alamouti_decoder_cc()
            sink = blocks.vector_sink_c()
            self.tb.connect(src, alamouti, sink)
            # Run flowgraph.
            self.tb.run()

            # Check if the expected result equals the actual result.
            self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

    ''' 
    5 tests validating the correct output of the decoder with random input data and random vector length.
    '''
    def test_002_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5
        num_tags = 4

        for i in range(repetitions):
            # Generate random input data with a random vector length.
            vlen = np.random.randint(2, 17)
            data = np.random.randn(vlen, data_length) + 1j * np.random.randn(vlen, data_length)
            # Generate random tag positions.
            tag_pos = np.random.randint(low=0, high=data_length/2, size=num_tags)*2
            tag_pos = np.sort(tag_pos)
            # Calculate expected result.
            tags, expected_result = self.dice_csi_tags(data,
                                                       num_tags,
                                                       tag_pos,
                                                       vlen)

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=np.reshape(data.T, vlen*data_length),
                                         repeat=False,
                                         tags=tags,
                                         vlen=vlen)

            alamouti = digital.alamouti_decoder_cc(vlen=vlen)
            sink = blocks.vector_sink_c()
            self.tb.connect(src, alamouti, sink)
            # Run flowgraph.
            self.tb.run()

            # Check if the expected result equals the actual result.
            self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

if __name__ == '__main__':
    gr_unittest.run(qa_alamouti_decoder_cc, "qa_alamouti_decoder_cc.xml")
