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

from matplotlib.pyplot import *

M_SQRT_2 = 1.0/np.sqrt(2)

class qa_diff_stbc_loopback (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    ''' 
    5 loopback encoder-decoder tests with random input, random modulation order,
    random channel matrix and random vector length.
    '''
    def test_001_t (self):
        # Define test params.
        data_length = 20
        repetitions = 5

        for n in range(repetitions):
            vlen = np.random.randint(1, 9)
            modulation_order = np.random.randint(1, 4)
            phase_shift = 2.0 * np.pi * np.random.randn()
            # Generate random input data.
            data = M_SQRT_2 * np.exp(1j* (2.0*np.pi*np.random.randint(0, 2**modulation_order, size=[data_length*vlen])/(2.0**modulation_order) + phase_shift))

            # Randomly generate normalized channel matrix.
            channel_gain_dist = np.random.rand()
            channel_matrix = np.array([channel_gain_dist*np.exp(2j*np.pi*np.random.rand()),
                                       np.sqrt(1-np.square(np.abs(channel_gain_dist)))*np.exp(2j*np.pi*np.random.rand())])

            # Set a stream tag to the beginning of the stream.
            tags = [(gr.tag_utils.python_to_tag((0,
                                                 pmt.string_to_symbol("packet_length"),
                                                 pmt.from_long(data_length),
                                                 pmt.from_long(0)))),

                    ]

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data, tags=tags)
            diff_stbc_encoder = digital.diff_stbc_encoder_cc(phase_shift, vlen)
            # Simulate channel with matrix multiplication.
            channel = blocks.multiply_matrix_cc_make([channel_matrix])
            v2s = blocks.stream_to_vector(gr.sizeof_gr_complex, vlen)
            diff_stbc_decoder = digital.diff_stbc_decoder_cc(phase_shift, vlen)
            sink = blocks.vector_sink_c()
            encoder_sink1 = blocks.vector_sink_c()
            encoder_sink2 = blocks.vector_sink_c()
            self.tb.connect(src, diff_stbc_encoder, channel, v2s, diff_stbc_decoder, sink)
            self.tb.connect((diff_stbc_encoder, 1), (channel, 1))
            self.tb.connect((diff_stbc_encoder, 0), encoder_sink1)
            self.tb.connect((diff_stbc_encoder, 1), encoder_sink2)
            # Run flowgraph.
            self.tb.run()

            self.assertComplexTuplesAlmostEqual(data, sink.data()[2*vlen:], 4)

if __name__ == '__main__':
    gr_unittest.run(qa_diff_stbc_loopback, "qa_diff_stbc_loopback.xml")
