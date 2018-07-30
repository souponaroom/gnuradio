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


class qa_ofdm_correct_carrier_freq_offset_vcvc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    '''
    Test validating the correct output data of the frequency corrector with
    random input data and 5 frames with random carrier frequency offset.'''
    def test_001_t (self):
        # Define test params.
        fft_len = 16
        frame_len = 4
        num_frames = 5
        tag_key = "carrier_freq_offset"
        # Produce random data.
        data = np.random.randn(fft_len*frame_len*num_frames) + 1j*np.random.randn(fft_len*frame_len*num_frames)
        # Random carrier freq offset.
        freq_off = np.random.randint(-fft_len, fft_len, num_frames)
        # Add tags at the start of each frame.
        tags = []
        for i in range(0, num_frames):
            tags.append(gr.tag_utils.python_to_tag((i*frame_len,
                                                    pmt.string_to_symbol(tag_key),
                                                    pmt.from_long(freq_off[i]),
                                                    pmt.from_long(0))))
        # Make gr blocks.
        src = blocks.vector_source_c(data=data, vlen=fft_len, tags=tags)
        carrier_correction = digital.ofdm_correct_carrier_freq_offset_vcvc(fft_len=fft_len,
                                                                           cp_len=0,
                                                                           carrier_freq_offset_key=tag_key)
        sink = blocks.vector_sink_c(vlen=fft_len)
        self.tb.connect(src, carrier_correction, sink)
        # Run flowgraph.
        self.tb.run ()
        # Calculate expected result
        expected_result = np.zeros(fft_len*frame_len*num_frames, dtype=complex)
        for i in range(0, num_frames):
            if freq_off[i] < 0:
                expected_result[i*fft_len*frame_len-freq_off[i]:(i+1)*fft_len*frame_len] = data[i*fft_len*frame_len:(i+1)*fft_len*frame_len+freq_off[i]]
            else:
                expected_result[i*fft_len*frame_len:(i+1)*fft_len*frame_len-freq_off[i]] = data[i*fft_len*frame_len+freq_off[i]:(i+1)*fft_len*frame_len]

        # Compare the expected result with the actual result.
        self.assertComplexTuplesAlmostEqual(expected_result, sink.data(), 2)

if __name__ == '__main__':
    gr_unittest.run(qa_ofdm_correct_carrier_freq_offset_vcvc, "qa_ofdm_correct_carrier_freq_offset_vcvc.xml")
