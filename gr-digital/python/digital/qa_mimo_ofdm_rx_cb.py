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
from gnuradio import analog
from gnuradio import fft
import digital_swig as digital
from gnuradio.digital.ofdm_txrx import ofdm_tx
from gnuradio.digital.mimo_ofdm_rx_cb import mimo_ofdm_rx_cb
import numpy as np

class qa_mimo_ofdm_rx_cb (gr_unittest.TestCase):
    '''
    4 basic loopback tests of the whole MIMO-OFDM transceiver with a
    random channel, no noise and different kinds of frequency offsets.
        '''
    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def simulate_loopback(self, data,
                          m, n, mimo_technique,
                          packet_len, packet_len_tag_key,
                          fft_len, cp_len,
                          f_off_rel):
        channel_matrix = (np.random.randn(n, m) + 1j * np.random.randn(n, m))

        # Source.
        src = blocks.vector_source_b(data, True, 1, ())
        s2tagged_stream = blocks.stream_to_tagged_stream(gr.sizeof_char, 1,
                                                         packet_len,
                                                         packet_len_tag_key)
        # MIMO-OFDM TX
        tx = ofdm_tx(
            fft_len=fft_len, cp_len=cp_len,
            packet_length_tag_key=packet_len_tag_key,
            bps_header=1,
            bps_payload=1,
            rolloff=0,
            debug_log=False,
            scramble_bits=False,
            m=m, mimo_technique=mimo_technique
        )
        # Static channel simulation.
        static_channel = blocks.multiply_matrix_cc(channel_matrix)
        # Apply frequency offset.
        const = analog.sig_source_c(fft_len, analog.GR_COS_WAVE, f_off_rel, 1.0)

        # MIMO-OFDM RX
        rx = mimo_ofdm_rx_cb(
            m=m, n=n,
            mimo_technique=mimo_technique,
            fft_len=fft_len,
            cp_len=cp_len,
            packet_length_tag_key=packet_len_tag_key,
            bps_header=1,
            bps_payload=1
        )
        sink = blocks.vector_sink_b()

        # Connect everything.
        self.tb.connect(src, blocks.head(gr.sizeof_char, packet_len * 100), s2tagged_stream, tx)
        for i in range(0, m):
            self.tb.connect((tx, i), (static_channel, i))
        mult = []
        for j in range(0, n):
            mult.append(blocks.multiply_cc())
            self.tb.connect((static_channel, j), mult[j], (rx, j))
            self.tb.connect(const, (mult[j], 1))
        self.tb.connect(rx, blocks.head(gr.sizeof_char, packet_len), sink)
        self.tb.run()
        
        return sink.data()


    def test_001_basic_t (self):
        """
        Basic test.
        """
        # Define test parameters.
        f_off_rel = 0.0
        packet_len = 8
        fft_len = 64
        cp_len = fft_len/4
        n = 2
        m = 2
        mimo_technique = "vblast"
        packet_len_tag_key = "packet_length"

        data = np.random.randint(0, 256, packet_len*6)
        result = self.simulate_loopback(data, m, n, mimo_technique,
                                        packet_len, packet_len_tag_key,
                                        fft_len, cp_len,
                                        f_off_rel)

        self.assertComplexTuplesAlmostEqual(data[:packet_len], result, 2)

    def test_002_fract_carr_freq_off_t (self):
        """
        Test with fractional carrier frequency offset.
        """
        # Define test parameters.
        f_off_rel = 0.5
        packet_len = 8
        fft_len = 64
        cp_len = fft_len/4
        n = 2
        m = 2
        mimo_technique = "vblast"
        packet_len_tag_key = "packet_length"

        data = range(packet_len)
        result = self.simulate_loopback(data, m, n, mimo_technique,
                                        packet_len, packet_len_tag_key,
                                        fft_len, cp_len,
                                        f_off_rel)

        self.assertComplexTuplesAlmostEqual(data, result, 2)


    def test_003_int_carr_freq_off_t (self):
        """
        Test with integer carrier frequency offset.
        """
        # Define test parameters.
        f_off_rel = 4.0
        packet_len = 8
        fft_len = 64
        cp_len = fft_len/4
        n = 2
        m = 2
        mimo_technique = "vblast"
        packet_len_tag_key = "packet_length"

        data = range(packet_len)
        result = self.simulate_loopback(data, m, n, mimo_technique,
                                        packet_len, packet_len_tag_key,
                                        fft_len, cp_len,
                                        f_off_rel)
        self.assertComplexTuplesAlmostEqual(data, result, 2)

    def test_004_arbitrary_freq_off_t (self):
        """
        Test with arbitrary frequency offset (int + fractional offset).
        """
        # Define test parameters.
        f_off_rel = 5.3
        packet_len = 8
        fft_len = 64
        cp_len = fft_len/4
        n = 2
        m = 2
        mimo_technique = "vblast"
        packet_len_tag_key = "packet_length"

        data = range(packet_len)
        result = self.simulate_loopback(data, m, n, mimo_technique,
                                        packet_len, packet_len_tag_key,
                                        fft_len, cp_len,
                                        f_off_rel)
        self.assertComplexTuplesAlmostEqual(data, result, 2)


if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_rx_cb, "qa_mimo_ofdm_rx_cb.xml")
