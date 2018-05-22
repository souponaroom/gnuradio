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
    def dice_csi_tags(self, input, num_tags, tag_pos):
        tags = []
        # Calculate initial behaviour before first tag.
        output = np.empty(shape=[len(input)], dtype=complex)
        output[::2] = (input[::2] + np.conj(input[1::2]))/2.0
        output[1::2] = (input[::2] - np.conj(input[1::2]))/2.0

        # Iterate over tags and update the calculated output according to the diced CSI.
        for i in range(0, num_tags):
            # Randomly generate CSI for one symbol.
            csi = (np.random.randn(2) + 1j * np.random.randn(2))
            # Assign the CSI vector to a PMT vector.
            csi_pmt = pmt.make_c32vector(2, csi[0])
            for k, channel in enumerate(csi):
                pmt.c32vector_set(v=csi_pmt, k=k, x=channel)
            # Append stream tags with CSI to data stream.
            tags.append(gr.tag_utils.python_to_tag((tag_pos[i],
                                                    pmt.string_to_symbol("csi"),
                                                    csi_pmt,
                                                    pmt.from_long(0))))
            # Calculate the expected result.
            total_branch_energy = np.sum(np.square(np.abs(csi)))
            if tag_pos[i]%2 == 0:
                output[tag_pos[i]  ::2] = (np.conj(csi[0])*input[tag_pos[i]::2] +
                                           csi[1]*np.conj(input[tag_pos[i]+1::2]))/total_branch_energy
                output[tag_pos[i]+1::2] = (np.conj(csi[1])*input[tag_pos[i]::2] -
                                           csi[0]*np.conj(input[tag_pos[i]+1::2]))/total_branch_energy
            else:
                output[tag_pos[i]+1::2] = (np.conj(csi[0])*input[tag_pos[i]+1::2] +
                                           csi[1]*np.conj(input[tag_pos[i]+2::2]))/total_branch_energy
                output[tag_pos[i]+2::2] = (np.conj(csi[1])*input[tag_pos[i]+1::2] -
                                           csi[0]*np.conj(input[tag_pos[i]+2::2]))/total_branch_energy

        return tags, output


    # 5 tests validating the correct output of the decoder with random input data.
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5
        num_tags = 4

        for i in range(repetitions):
            # Generate random input data.
            data = np.random.randn(data_length) + 1j * np.random.randn(data_length)
            # Generate random tag positions.
            tag_pos = np.random.randint(low=0, high=data_length/2, size=num_tags)*2
            tag_pos = np.sort(tag_pos)
            # Calculate expected result.
            tags, expected_result = self.dice_csi_tags(data,
                                                       num_tags,
                                                       tag_pos)

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data,
                                         repeat=False,
                                         tags=tags)
            alamouti = digital.alamouti_decoder_cc()
            sink = blocks.vector_sink_c()
            self.tb.connect(src, alamouti, sink)
            # Run flowgraph.
            self.tb.run()

            # Check if the expected result equals the actual result.
            self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)

if __name__ == '__main__':
    gr_unittest.run(qa_alamouti_decoder_cc, "qa_alamouti_decoder_cc.xml")
