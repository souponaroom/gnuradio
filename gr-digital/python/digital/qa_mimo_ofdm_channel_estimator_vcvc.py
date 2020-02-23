#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018,2019,2020 Moritz Luca Schmid, 
# Communications Engineering Lab (CEL) / Karlsruhe Institute of Technology (KIT). 
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
from gnuradio import fft
import digital_swig as digital
from gnuradio.digital.ofdm_txrx import ofdm_tx
from mimo import mimo_technique as mimo
import numpy as np
import pmt

class qa_mimo_ofdm_channel_estimator_vcvc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    """
    This unit test validates the correct channel estimation of a constant and flat channel
    with random coefficients.
    TODO enable N neq M
    """
    def test_static_channel_estimation_t (self):
        # Channel dimesions to test.
        dimensions = [1, 2, 4, 8]
        for dim in dimensions:
            self.estimate_NxN_channel(dim)

    def estimate_NxN_channel (self, channel_dimension):
        # Define test params.
        N = channel_dimension
        packet_len = N # must be a multiple of N (so that min one whole pilot sequence fits in each packet) 
        len_tag_key = 'packet_length'
        fft_len = 4
        pilot_carriers = [-fft_len/2, fft_len/2-1] # use the outer carriers as pilots
        occupied_carriers = range(-fft_len/2+1, fft_len/2-1) # all other carriers are occupied with payload
        channel_matrix = (np.random.randn(N, N) + 1j * np.random.randn(N, N))
        self.walsh_sequences = np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, -1, 1, -1, 1, -1, 1, -1],
            [1, 1, -1, -1, 1, 1, -1, -1],
            [1, -1, -1, 1, 1, -1, -1, 1],
            [1, 1, 1, 1, -1, -1, -1, -1],
            [1, -1, 1, -1, -1, 1, -1, 1],
            [1, 1, -1, -1, -1, -1, 1, 1],
            [1, -1, -1, 1, -1, 1, 1, -1]
        ])

        # prepare data for antenna 1 & 2
        data = np.repeat(self.walsh_sequences[:N,:N],fft_len, axis=1)

        # Static channel.
        static_channel = blocks.multiply_matrix_cc(channel_matrix)

        # Define TX blocks.
        src = [] # data sources
        ofdm_mod = [] # ifft for each TX branch
        v2s = [] # vector to stream (to feed it to the channel matrix)
        for i in range(N):
            src.append( blocks.vector_source_c(data[i], False, fft_len, ()))
            ofdm_mod.append( fft.fft_vcc(fft_len, False, (), True))
            v2s.append( blocks.vector_to_stream(gr.sizeof_gr_complex, fft_len))
        for i in range(N):
            # Connect.
            self.tb.connect(src[i], ofdm_mod[i], v2s[i], (static_channel, i))

        # Channel estimator (the target of this test).
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(m=N, n=N,
                                                               fft_len=fft_len,
                                                               pilot_symbols=self.walsh_sequences[:N, :N],
                                                               pilot_carriers=pilot_carriers,
                                                               occupied_carriers=occupied_carriers,
                                                               csi_key="csi", start_key="start")

        # Define RX blocks.
        add_tag = [] # add packet tag for each RX branch
        ofdm_demod = [] # FFT as OFDM demod for each RX branch
        s2v = [] # stream to vector (to feed to estimator)
        norm = [] # normalize each RX branch
        sink = [] # data sink for each RX branch (after estimation)
        for i in range(N):
            add_tag.append( blocks.stream_to_tagged_stream(gr.sizeof_gr_complex, 
                                                           fft_len,
                                                           packet_len,
                                                           len_tag_key))
            ofdm_demod.append( fft.fft_vcc(fft_len, True, (), True))
            s2v.append( blocks.stream_to_vector(gr.sizeof_gr_complex, fft_len))
            norm.append( blocks.multiply_const_cc(1./fft_len, fft_len))
            sink.append(blocks.vector_sink_c(vlen=fft_len-len(pilot_carriers)))
            # Connect.
            self.tb.connect((static_channel,i), s2v[i], ofdm_demod[i], norm[i], 
                            add_tag[i], (channel_est, i), sink[i])
            
        
        # Run test.
        self.tb.run ()

        # Read channel estimation (out of tags) and compare with actual channel state.
        csi = np.empty(shape=[len(occupied_carriers), N, N], dtype=complex)
        for k in range(0, len(occupied_carriers)):
            for n in range(0, N):
                for m in range(0, N):
                    csi[k][n][m] = pmt.c32vector_ref(pmt.vector_ref(pmt.vector_ref(sink[0].tags()[0].value, k), n), m)

        for c in range(0, len(occupied_carriers)):
            for n in range(0, N):
                self.assertComplexTuplesAlmostEqual(channel_matrix[n], csi[c][n], 3)

if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_channel_estimator_vcvc, "qa_mimo_ofdm_channel_estimator_vcvc.xml")
