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
import digital_swig as digital
from gnuradio.digital.ofdm_txrx import ofdm_tx
from gnuradio.digital.mimo_ofdm_rx_cb import mimo_ofdm_rx_cb
from mimo import mimo_technique as mimo
import numpy as np
import pmt

'''
This qa test serves 2 functions:
1.) Generic loopback test for MIMO-OFDM transmission and reception.
    -flowgraph: random_data - mimo_ofdm_tx - static_channel - mimo_ofdm_rx
    -generic parameters: -dimension

2.) Tutorial to 'How to build up a basic MIMO-OFDM transceiver' with 
    comments and practical explanations to each step.
'''
class qa_mimo_ofdm_loopback (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_mimo_ofdm_loopback_t (self):
        # Channel dimesions to test.
        dimensions = [1, 2, 4, 8]
        for dim in dimensions:
            self.parametrized_loopback_test(dim)

    def parametrized_loopback_test (self, channel_dimension):
        print("dimension")
        print(channel_dimension)
        # Define test params.
        N = channel_dimension
        mimo_rate = N # like a code rate the number of parallel streams per time slot
        fft_len = 64 
        cp_len = fft_len/8
        # Number of OFDM symbols per MIMO branch per packet.
        # (including data for header, payload and CRC)
        packet_len = 3+N 
        pilot_carriers = [[-12, -4, 4, 12]] # We chose 4 equally spaced pilot carriers.
        # All other carriers are occupied with payload.
        occupied_carriers = [[x for x in range(-fft_len/4-2, fft_len/4+2) 
                                if x not in pilot_carriers[0]]] 
        channel_matrix = (np.random.randn(N, N) + 1j * np.random.randn(N, N))

        # create random bytes and set up source blocks
        payload_bits_per_packet = len(occupied_carriers[0])*(mimo_rate*packet_len-1)-32 #
        assert payload_bits_per_packet > 0, "The number of payload bits per packet (%d) must be >0!" % payload_bits_per_packet
        assert payload_bits_per_packet % 8 == 0, "The number of bits per packet (%d) must be a multiple of 8!" % payload_bits_per_packet
        print("payload bits per packet")
        print(payload_bits_per_packet)
        payload_bytes_per_packet = payload_bits_per_packet/8
        data = np.random.randint(256, size=1*  payload_bytes_per_packet) # 3 packet
        src = blocks.vector_source_b(data, False, 1, ())
        add_tag = blocks.stream_to_tagged_stream(gr.sizeof_char, 1,
                                                 payload_bytes_per_packet, "packet_length")

        # MIMO-OFDM TX
        # as a pilot sequence we use Walsh sequences (default of these hier blocks)
        tx = ofdm_tx(m=N,
                     fft_len=fft_len,
                     cp_len=cp_len,
                     occupied_carriers=occupied_carriers,
                     pilot_carriers=pilot_carriers,
                     bps_header=1, bps_payload=1,
                     sync_word1=None, sync_word2=None,
                     frame_len_tag_key="frame_length",
                     packet_len_tag_key="packet_length",
                     packet_num_tag_key="packet_num",
                     rolloff=0, debug_log=False,
                     scramble_bits=False,
                     mimo_technique=mimo.VBLAST_ZF)

        # Static channel.
        static_channel = blocks.multiply_matrix_cc(channel_matrix)


        # MIMO-OFDM RX
        rx = mimo_ofdm_rx_cb(m=N, n=N,
                             mimo=mimo.VBLAST_ZF,
                             fft_len=fft_len, cp_len=cp_len,
                             frame_length_tag_key="frame_length",
                             packet_length_tag_key="packet_length",
                             packet_num_tag_key="packet_num",
                             occupied_carriers=occupied_carriers,
                             pilot_carriers=pilot_carriers,
                             bps_header=1, bps_payload=1,
                             sync_word1=None, sync_word2=None,
                             scramble_bits=False, show_const=False)
        # Connect.
        self.tb.connect(src, add_tag, tx)

        tx_sinks = []
        for i in range(N):
            tx_sinks.append(blocks.vector_sink_c())
            self.tb.connect((tx,i), (static_channel,i))
            self.tb.connect((static_channel,i), tx_sinks[i])
        
        self.tb.run()

        rx_data = []
        rx_srcs = []
        for i in range(N):
            assert len(tx_sinks[i].data()) ==1* (cp_len+fft_len)*(packet_len+2) # 2 packets
            rx_data.append(np.append(np.append( np.random.randn((fft_len+cp_len)*4), tx_sinks[i].data()), np.random.randn((fft_len+cp_len)*10)))
            #rx_data.append(tx_sinks[i].data())
            rx_srcs.append(blocks.vector_source_c(np.append(rx_data[i],rx_data[i]), False, 1, ()))
            self.tb.connect(rx_srcs[i], (rx, i))


        # Sink.
        sink = blocks.vector_sink_b()
        self.tb.connect(rx, sink)
            
        # Run test.
        self.tb.run()

        # Compare input and output data.
     
    #print(data)
    #    print(sink.data())
        self.assertComplexTuplesAlmostEqual(np.append(data,data), sink.data(), 2)


if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_loopback, "qa_mimo_ofdm_loopback.xml")
