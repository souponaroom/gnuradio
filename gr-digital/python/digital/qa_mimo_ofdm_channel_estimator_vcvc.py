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

import numpy as np
from scipy.linalg import hadamard
import pmt

"""
This qa unit test validates the block mimo_ofdm_channel_estimator_vcvc.
The block is used for the channel estimation and a time delay correction of a
MIMO-OFDM receiver. 

In detail, the following units are tested:
-Estimation of a 
    - constant and static channel with random coefficients.
    - constant but time variant channel with random coefficients.
-Estimation and correction of a cyclic time shift, 
 which can happen at the MIMO-OFDM time synchronization.

Setup: 
src - ofdm_mod(fft) - static_channel - ofdm_demod(ifft) - estimator(to test) - sink

Notes:
-The packet length must be >= the MIMO dimension (to separate the channels).
-This qa test uses 3 packets of different lenghts to verify correct behavior.
-The unit tests use the following MIMO dimensions: 1, 2, 4.
"""
class qa_mimo_ofdm_channel_estimator_vcvc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    # Unit test for constant and static channel estimation.
    def test_constant_static_channel_estimation_t (self):
        for N in [1,2,4]:
            channel, csi, tag_pos, delay = self.build_and_run_flowgraph(N, static=True)
            [[np.testing.assert_almost_equal(np.array(channel), np.array(carrier_est), 3) 
             for carrier_est in symbol] for symbol in csi]

    # Unit test for constant and time variant channel estimation.
    def test_constant_time_variant_channel_estimation_t (self):
        for N in [1,2,4]:
            channel, csi, tag_pos, delay = self.build_and_run_flowgraph(N, static=False)
            [[np.testing.assert_almost_equal(
                    np.array(channel[np.nonzero(tag_pos <= i)[-1][-1]]), 
                    np.array(carrier_est), 3) 
              for carrier_est in symbol] 
             for i,symbol in enumerate(csi)]

    # Unit test for estimation and correction of a cyclic time shift.
    def test_cyclic_time_shift_estimation_t (self):
        # Test 3 arbitrary combinations of MIMO dimension and delay.
        for (N, time_shift) in [(2,0),(1,1),(4,3)]: 
            channel, csi, tag_pos, delay = self.build_and_run_flowgraph(N, 
                                                                        static=True, 
                                                                        delay=time_shift)
            self.assertEqual(time_shift, delay)

    def build_and_run_flowgraph (self, N, static, delay=0):
        # Define flowgraph parameters.
        fft_len = 4
        # Use the outer carriers as pilots.
        pilot_carriers = [-fft_len/2, fft_len/2-1] 
        # All other carriers (2) are occupied with payload.
        occupied_carriers = range(-fft_len/2+1, fft_len/2-1) 

        # Generate pilots by using Walsh sequences. 
        # The payload data is irrelevant here and thus chosen to be 
        # equal to the pilots to keep it simple.
        data = np.repeat(hadamard(N),fft_len, axis=1)

        # Test different packet lengths N/2, N and 7N/2 (must be >= N to separate MIMO channels). 
        tag_pos = np.array([0, N+N/2, 2*N+N/2, 6*N])

        if static:
            # Generate static random channel matrix here.
            channel_matrix = (np.random.randn(N, N) + 1j * np.random.randn(N, N))
        else:
            channel_matrix = (np.random.randn(len(tag_pos)-1, N, N) + 1j * 
                              np.random.randn(len(tag_pos)-1, N, N))

        channel_data = np.empty((N, fft_len*tag_pos[-1]), dtype=complex)

        tags = []
        # Apply channel on each packet (with varying channel matrix if static=False)
        for c in range(len(tag_pos)-1):
            tags.append(gr.tag_utils.python_to_tag((tag_pos[c]*fft_len,
                                                    pmt.string_to_symbol("start"),
                                                    pmt.from_long(0))))
            channel = blocks.multiply_matrix_cc(channel_matrix if static else channel_matrix[c])

            src = [blocks.vector_source_c(np.tile(data[i],6)[:fft_len*(tag_pos[c+1]-tag_pos[c])], 
                                          False, fft_len) for i in range(N)]
            ofdm_mod = [fft.fft_vcc(fft_len, False, (), True) for i in range(N)]
            v2s = [blocks.vector_to_stream(gr.sizeof_gr_complex, fft_len) for i in range(N)]
            sink = [blocks.vector_sink_c() for i in range(N)]

            # Connect TX blocks with channel.
            for i in range(N): self.tb.connect(src[i], ofdm_mod[i], v2s[i], (channel, i), sink[i])

            # Run flowgraph and apply channel on each packet.
            self.tb.run ()
            # Save results of this packet..
            channel_data[:,fft_len*tag_pos[c]:fft_len*tag_pos[c+1]] = np.array([x.data() for x in sink])


        # Apply cyclic time shift to simulate a non perfect time synchonization.
        channel_data = channel_data[:,np.roll(np.arange(fft_len*tag_pos[-1]).
                       reshape((-1,fft_len)),delay, axis=1).flatten()]
        
        # Channel estimator (the target of this test).
        channel_est = digital.mimo_ofdm_channel_estimator_vcvc(
                        m=N, n=N, fft_len=fft_len, pilot_symbols=hadamard(N),
                        pilot_carriers=pilot_carriers, 
                        occupied_carriers=occupied_carriers,
                        csi_key="csi", start_key="start")

        # Define RX blocks.
        src2 = [blocks.vector_source_c(channel_data[i], False, 1, tags=tags) for i in range(N)]
        ofdm_demod = [fft.fft_vcc(fft_len, True, (), True) for i in range(N)]
        s2v = [blocks.stream_to_vector(gr.sizeof_gr_complex, fft_len) for i in range(N)]
        norm = [blocks.multiply_const_cc(1./fft_len, fft_len) for i in range(N)]
        sink = [blocks.vector_sink_c(vlen=fft_len-len(pilot_carriers)) for i in range(N)]

        # Connect channel estimator to flowgraph.
        for i in range(N): self.tb.connect(src2[i], s2v[i], ofdm_demod[i], norm[i], 
                            (channel_est, i), sink[i])
            
        # Run test.
        self.tb.run()

        csi_tags = [x for x in sink[0].tags() if pmt.symbol_to_string(x.key)=="csi"]

        csi = []
        for i in range(tag_pos[-1]):
            # Read channel estimation (out of tags) and compare with actual channel state.
            csi.append([[[pmt.c32vector_ref(pmt.vector_ref(pmt.vector_ref(csi_tags[i].value, k), n), m) 
                     for m in range(N)] 
                    for n in range(N)] 
                   for k in range(len(occupied_carriers))])

        return channel_matrix, csi, tag_pos, channel_est.get_time_delay()


if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_channel_estimator_vcvc, "qa_mimo_ofdm_channel_estimator_vcvc.xml")
