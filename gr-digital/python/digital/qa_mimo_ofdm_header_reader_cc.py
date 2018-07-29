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
from gnuradio import fft
import digital_swig as digital
from gnuradio.digital.ofdm_txrx import ofdm_tx
import numpy as np
import pmt

class qa_mimo_ofdm_header_reader_cc (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def _get_constellation(self, bps):
        """ Returns a modulator block for a given number of bits per symbol """
        constellation = {
            1: digital.constellation_bpsk(),
            2: digital.constellation_qpsk(),
            3: digital.constellation_8psk()
        }
        try:
            return constellation[bps]
        except KeyError:
            print 'Modulation not supported.'
            exit(1)

    ''' Generate header and modulate payload data.'''
    def generate_tx_data(self, data, packet_lengths,
                         occupied_carriers,
                         bps_header, bps_payload,
                         packet_len_tag_key,
                         header_constellation, payload_constellation):
        # Append stream tags to data.
        tags = []
        offset = 0
        for i in range(0, packet_lengths.size):
            tags.append(gr.tag_utils.python_to_tag((offset,
                                                    pmt.string_to_symbol(packet_len_tag_key),
                                                    pmt.from_long(packet_lengths[i]))))
            offset += packet_lengths[i]
        payload_src = blocks.vector_source_b(data=data, tags=tags)

        # Generate header.
        formatter_object_tx = digital.packet_header_ofdm(
            occupied_carriers=[occupied_carriers], n_syms=1,
            bits_per_header_sym=bps_header,
            bits_per_payload_sym=bps_payload,
            scramble_header=False)
        header_gen = digital.packet_headergenerator_bb(formatter_object_tx.base(),
                                                       packet_len_tag_key)

        header_mod = digital.chunks_to_symbols_bc(header_constellation.points())
        header_payload_mux = blocks.tagged_stream_mux(
            itemsize=gr.sizeof_gr_complex * 1,
            lengthtagname=packet_len_tag_key,
            tag_preserve_head_pos=1)  # Head tags on the payload stream stay on the head
        self.tb.connect(header_gen, header_mod, (header_payload_mux, 0))

        # Payload modulation
        payload_unpack = blocks.repack_bits_bb(
            8,  # Unpack 8 bits per byte
            bps_payload,
            packet_len_tag_key)

        payload_mod = digital.chunks_to_symbols_bc(payload_constellation.points())

        self.tb.connect(payload_src, payload_unpack, payload_mod, (header_payload_mux, 1))
        self.tb.connect(payload_src, header_gen)
        sink = blocks.vector_sink_c()
        self.tb.connect(header_payload_mux, sink)
        self.tb.run()
        return sink.data()

    '''Accepts tagged stream data and reads/parses the header and demodulates the payload.'''
    def demod(self, rx_data, rx_tags,
              occupied_carriers,
              bps_header, bps_payload,
              header_constellation, payload_constellation
              ):

        rx_src = blocks.vector_source_c(data=rx_data,
                                        tags=rx_tags)
        # Header reader.
        header_formatter_rx = digital.packet_header_ofdm(
            [occupied_carriers], 1,
            "packet_length",
            "frame_length",
            "packet_num",
            bps_header, bps_payload)
        header_reader = digital.mimo_ofdm_header_reader_cc(header_constellation.base(),
                                                           header_formatter_rx.formatter())
        # Payload demodulator.
        payload_demod = digital.constellation_decoder_cb(payload_constellation.base())
        payload_pack = blocks.repack_bits_bb(bps_payload, 8, "", True)
        sink_rx = blocks.vector_sink_b()
        self.tb.connect(rx_src, header_reader, payload_demod, payload_pack, sink_rx)
        self.tb.run()
        return sink_rx.data()

    '''Test with 2 packets which follow each other directly.'''
    def test_001_t (self):
        # Define test params.
        packet_lengths = np.array([4, 4])
        data_length = np.sum(packet_lengths)
        occupied_carriers = range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22, 27)
        header_len = len(occupied_carriers)
        bps_header = 1
        bps_payload = 1
        packet_len_tag_key = "packet_length"

        # Constellations.
        header_constellation = self._get_constellation(bps_header)
        payload_constellation = self._get_constellation(bps_payload)

        # Payload data.
        data = np.random.randint(0, 100, data_length)

        # Generate tx data.
        tx_data = self.generate_tx_data(data, packet_lengths,
                                        occupied_carriers,
                                        bps_header, bps_payload,
                                        packet_len_tag_key,
                                        header_constellation, payload_constellation)

        # Manipulate data and tags
        rx_data = tx_data

        # Read data from sink and add 'start' tags.
        rx_tags = []
        offset = 0
        for i in range(0, packet_lengths.size):
            rx_tags.append(gr.tag_utils.python_to_tag((offset,
                                                       pmt.string_to_symbol("start"),
                                                       pmt.from_long(0))))
            offset += packet_lengths[i] * 8 / bps_payload + header_len

        # Header and payload reader.
        rx_data = self.demod(rx_data, rx_tags,
                             occupied_carriers,
                             bps_header, bps_payload,
                             header_constellation, payload_constellation)

        # Check if the expected result equals the actual result.
        self.assertComplexTuplesAlmostEqual(data, rx_data, 2)

if __name__ == '__main__':
    gr_unittest.run(qa_mimo_ofdm_header_reader_cc, "qa_mimo_ofdm_header_reader_cc.xml")
