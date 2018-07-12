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
from gnuradio import digital
import numpy as np
import pmt

class qa_mimo_encoder_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def run_test(self, M, mimo_technique):
        payload_length = 6
        payload = np.random.randn(payload_length) + 1j * np.random.randn(payload_length)
        # Produce tags.
        tags = [gr.tag_utils.python_to_tag((0,
                                            pmt.string_to_symbol("length"),
                                            pmt.from_long(payload_length),
                                            pmt.from_long(0)))]
        # Build GNU Radio blocks.
        training_length = 4
        training_sequence = np.random.randn(M, training_length)

        src = blocks.vector_source_c(data=np.append(training_sequence[0], payload),
                                      repeat=False,
                                      tags=tags)
        encoder = digital.mimo_encoder_cc(M=M,
                                          mimo_technique=mimo_technique,
                                          length_tag_name='length',
                                          training_sequence=training_sequence)
        self.tb.connect(src, encoder)
        sink = []
        for m in range(0, M):
            sink.append(blocks.vector_sink_c_make())
            self.tb.connect((encoder, m), sink[m])
        # Run flowgraph.
        self.tb.run()
        # Check if the training sequences were inserted correctly at each outport.
        for m in range(0, M):
            self.assertComplexTuplesAlmostEqual(training_sequence[m], sink[m].data()[0:len(training_sequence)], 2)


    def test_001_t (self):
        # Test parameters.
        M = 2
        repetitions = 2
        for i in range(0, repetitions):
            self.run_test(M, 'alamouti')



if __name__ == '__main__':
    gr_unittest.run(qa_mimo_encoder_cc, "qa_mimo_encoder_cc.xml")
