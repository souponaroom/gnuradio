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

class qa_alamouti_loopback (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # 5 loopback encoder-decoder tests with random input data and random CSI.
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5
        num_tags = 4
        result = np.empty(shape=[data_length], dtype=complex)

        for n in range(repetitions):
            # Generate random input data.
            data = np.random.randn(data_length) + 1j * np.random.randn(data_length)
            # Generate random tag positions.
            tag_pos = np.random.randint(low=0, high=data_length/2, size=num_tags)*2.0
            # Add tag pos 0 for initial channel state.
            tag_pos = np.append(tag_pos, 0)
            tag_pos = np.sort(tag_pos)

            # Iterate over tags.
            for i in range(0, num_tags):
                # Randomly generate CSI for one symbol.
                csi = (np.random.randn(2) + 1j * np.random.randn(2))
                # Assign the CSI vector to a PMT vector.
                csi_pmt = pmt.make_c32vector(2, csi[0])
                for k, channel in enumerate(csi):
                    pmt.c32vector_set(v=csi_pmt, k=k, x=channel)
                # Append stream tags with CSI to data stream.
                tags = [(gr.tag_utils.python_to_tag((0,
                                                     pmt.string_to_symbol("csi"),
                                                     csi_pmt,
                                                     pmt.from_long(0))))]

                # Build up the test flowgraph.
                src = blocks.vector_source_c(data=data[tag_pos[i]::],
                                             repeat=False,
                                             tags=tags)
                alamouti_encoder = digital.alamouti_encoder_cc()
                # Simulate channel with matrix multiplication.
                channel = blocks.multiply_matrix_cc_make([csi])
                alamouti_decoder = digital.alamouti_decoder_cc()
                sink = blocks.vector_sink_c()
                self.tb.connect(src, alamouti_encoder, channel, alamouti_decoder, sink)
                self.tb.connect((alamouti_encoder, 1), (channel, 1))
                # Run flowgraph.
                self.tb.run()
                result[tag_pos[i]::] = sink.data()

            ''' 
            Check if the expected result (=the data itself, because 
            we do a loopback) equals the actual result. '''
            self.assertComplexTuplesAlmostEqual(data, result, 4)

if __name__ == '__main__':
    gr_unittest.run(qa_alamouti_loopback, "qa_alamouti_loopback.xml")
