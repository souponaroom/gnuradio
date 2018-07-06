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

M_SQRT_2 = 1.0/np.sqrt(2)

class qa_diff_stbc_decoder_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Produce stream tags and calculate expected result.
    def decode(self, basis, input, tag_pos):
        # Calculate the expected behaviour without the presence of tags.
        output = np.empty(shape=[len(input)], dtype=complex)
        data = np.append(basis, input)
        output[::2] = (data[2::2]*np.conj(data[0:len(input):2]) + np.conj(data[3::2])*data[1:len(input):2])*basis[0] - \
                      (data[2::2]*np.conj(data[1:len(input):2]) - np.conj(data[3::2])*data[0:len(input):2])*np.conj(basis[1])
        output[1::2] = (data[2::2]*np.conj(data[0:len(input):2]) + np.conj(data[3::2])*data[1:len(input):2])*basis[1] + \
                       (data[2::2]*np.conj(data[1:len(input):2]) - np.conj(data[3::2])*data[0:len(input):2])*np.conj(basis[0])

        # Iterate over tags and update the calculated output according to the random CSI.
        tags = []
        for i in range(0, len(tag_pos)):
            # Assign the CSI vector to a PMT vector.
            csi_pmt = pmt.from_long(True)
            # Append stream tags to data stream.
            tags.append(gr.tag_utils.python_to_tag((tag_pos[i],
                                                    pmt.string_to_symbol("start"),
                                                    csi_pmt,
                                                    pmt.from_long(0))))
        return tags, np.delete(output, np.append(tag_pos, tag_pos+1))

    ''' 5 test with random input data, random tag positions, random basis.'''
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5
        num_tags = 4

        for i in range(repetitions):
            data = np.random.randint(-1, 2, size=data_length) + 1j*np.random.randint(-1, 2, size=data_length)
            # Generate random tag positions.
            tag_pos = np.random.randint(low=0, high=data_length / 2, size=num_tags) * 2
            tag_pos = np.sort(tag_pos)

            phase_shift = 2.0 * np.pi * np.random.randn()
            basis = np.array([M_SQRT_2 * np.exp(1j * phase_shift), M_SQRT_2 * np.exp(1j * phase_shift)])

            tags, expected_result = self.decode(basis, data, tag_pos)

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data,
                                         repeat=False,
                                         tags=tags)
            stbc = digital.diff_stbc_decoder_cc(phase_shift)
            sink = blocks.vector_sink_c()
            self.tb.connect(src, stbc, sink)
            # Run flowgraph.
            self.tb.run()

            self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 4)


if __name__ == '__main__':
    gr_unittest.run(qa_diff_stbc_decoder_cc, "qa_diff_stbc_decoder_cc.xml")
