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

class qa_vblast_loopback (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def build_and_run_flowgraph(self, repetitions, data_length, num_inputs_min, num_inputs_max, equalizer_type, vlen=[1, 1]):
        for a in range(repetitions):
            num_inputs = np.random.randint(low=num_inputs_min, high=num_inputs_max+1)
            # Generate random input data.
            data = np.random.randn(data_length*num_inputs*vlen[a]) + 1j * np.random.randn(data_length*num_inputs*vlen[a])

            # Randomly generate CSI for one symbol.
            csi = (np.random.randn(vlen[a], num_inputs, num_inputs) + 1j * np.random.randn(vlen[a], num_inputs, num_inputs))

            # Assign the CSI vector to a PMT vector.
            csi_pmt = pmt.make_vector(vlen[a], pmt.make_vector(num_inputs, pmt.make_c32vector(num_inputs, 1.0)))
            for k, carrier in enumerate(csi):
                carrier_vector_pmt = pmt.make_vector(num_inputs, pmt.make_c32vector(num_inputs, csi[k][0][0]))
                for l, rx in enumerate(csi[k]):
                    line_vector_pmt = pmt.make_c32vector(num_inputs, csi[k][l][0])
                    for m, tx in enumerate(csi[k][l]):
                        pmt.c32vector_set(v=line_vector_pmt, k=m, x=csi[k][l][m])
                    pmt.vector_set(carrier_vector_pmt, l, line_vector_pmt)
                pmt.vector_set(csi_pmt, k, carrier_vector_pmt)

            # Append stream tags with CSI to data stream.

            tags = [(gr.tag_utils.python_to_tag((0,
                                                 pmt.string_to_symbol("csi"),
                                                 csi_pmt,
                                                 pmt.from_long(0))))]

            if equalizer_type == 'MMSE':
                # Add an SNR tag at the start of the stream for MMSE.
                tags.append(gr.tag_utils.python_to_tag((0,
                                                        pmt.string_to_symbol("snr"),
                                                        pmt.make_f32vector(num_inputs, 1e8),
                                                        pmt.from_long(0))))

            # Build up the test flowgraph.
            src = blocks.vector_source_c(data=data,
                                         repeat=False,
                                         tags=tags)
            vblast_encoder = digital.vblast_encoder_cc(num_inputs)
            demux = []
            channels = []
            for i in range(0, vlen[a]):
                for j in range(0, num_inputs):
                    demux.append(blocks.keep_m_in_n(gr.sizeof_gr_complex, 1, vlen[a], i))
                # Simulate channel with matrix multiplication.
                channels.append(blocks.multiply_matrix_cc_make(csi[i]))
            mux = []
            s2v = []
            for i in range(0, num_inputs):
                mux.append(blocks.stream_mux(gr.sizeof_gr_complex, [1]*vlen[a]))
                s2v.append(blocks.stream_to_vector(gr.sizeof_gr_complex, vlen[a]))

            vblast_decoder = digital.vblast_decoder_cc(num_inputs, equalizer_type, vlen[a])
            sink = blocks.vector_sink_c()
            self.tb.connect(src, vblast_encoder)
            for n in range(0, num_inputs):
                for i in range(0, vlen[a]):
                    self.tb.connect((vblast_encoder, n), demux[i*num_inputs+n], (channels[i], n))
                    self.tb.connect((channels[i], n), (mux[n], i))
                self.tb.connect(mux[n], s2v[n], (vblast_decoder, n))

            self.tb.connect(vblast_decoder, sink)
            # Run flowgraph.
            self.tb.run()

            ''' 
            Check if the expected result (=the data itself, because 
            we do a loopback) equals the actual result. '''
            self.assertComplexTuplesAlmostEqual(data, sink.data(), 2)

    '''
    2 tests validating the correct output of the loopback with random input data, ZF equalizer,
    1x1 MIMO scheme and vector lengths {1,[2,16]}. '''

    def test_001_t(self):
        self.build_and_run_flowgraph(repetitions=2,
                                     data_length=4,
                                     num_inputs_min=1,
                                     num_inputs_max=1,
                                     equalizer_type='ZF',
                                     vlen=[1, np.random.randint(2,17)])

    '''
    2 tests validating the correct output of the loopback with random input data, MMSE equalizer,
    1x1 MIMO scheme, vector lengths {1,[2,16]} and an extremely high SNR regime (1/snr -> 0). '''

    def test_002_t(self):
        self.build_and_run_flowgraph(repetitions=2,
                                     data_length=4,
                                     num_inputs_min=1,
                                     num_inputs_max=1,
                                     equalizer_type='MMSE',
                                     vlen=[1, np.random.randint(2, 17)])

    '''
    2 tests validating the correct output of the loopback with random input data, ZF equalizer,
    1x1 MIMO scheme and vector lengths {1,[2,16]}. '''

    def test_003_t(self):
        self.build_and_run_flowgraph(repetitions=2,
                                     data_length=4,
                                     num_inputs_min=2,
                                     num_inputs_max=2,
                                     equalizer_type='ZF',
                                     vlen=[1, np.random.randint(2, 17)])

    '''
    2 tests validating the correct output of the loopback with random input data, MMSE equalizer,
    1x1 MIMO scheme, vector lengths {1,[2,16]} and an extremely high SNR regime (1/snr -> 0). '''

    def test_004_t(self):
        self.build_and_run_flowgraph(repetitions=2,
                                     data_length=4,
                                     num_inputs_min=2,
                                     num_inputs_max=2,
                                     equalizer_type='MMSE',
                                     vlen=[1, np.random.randint(2, 17)])

    '''
    2 tests validating the correct output of the loopback with random input data, ZF equalizer,
    1x1 MIMO scheme and vector lengths {1,2}. '''

    def test_005_t(self):
        self.build_and_run_flowgraph(repetitions=2,
                                     data_length=4,
                                     num_inputs_min=3,
                                     num_inputs_max=16,
                                     equalizer_type='ZF',
                                     vlen=[1, 2])

    '''
    2 tests validating the correct output of the loopback with random input data, MMSE equalizer,
    1x1 MIMO scheme, vector lengths {1,2} and an extremely high SNR regime (1/snr -> 0). '''

    def test_006_t(self):
        self.build_and_run_flowgraph(repetitions=2,
                                     data_length=4,
                                     num_inputs_min=3,
                                     num_inputs_max=16,
                                     equalizer_type='MMSE',
                                     vlen=[1, 2])

if __name__ == '__main__':
    gr_unittest.run(qa_vblast_loopback, "qa_vblast_loopback.xml")
