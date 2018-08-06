/* -*- c++ -*- */
/* 
 * Copyright 2018 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/format.hpp>
#include <gnuradio/io_signature.h>
#include "mimo_ofdm_header_reader_cc_impl.h"

using namespace boost;

namespace gr {
  namespace digital {

    const pmt::pmt_t mimo_ofdm_header_reader_cc_impl::d_key = pmt::string_to_symbol("start");

    mimo_ofdm_header_reader_cc::sptr
    mimo_ofdm_header_reader_cc::make(constellation_sptr constellation,
                                     const gr::digital::packet_header_default::sptr &header_formatter)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_header_reader_cc_impl(constellation, header_formatter));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_header_reader_cc_impl::mimo_ofdm_header_reader_cc_impl(
            constellation_sptr constellation,
            const gr::digital::packet_header_default::sptr &header_formatter)
      : gr::block("mimo_ofdm_header_reader_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_constellation(constellation),
        d_dim(constellation->dimensionality()),
        d_header_formatter(header_formatter),
        d_header_length(header_formatter->header_len()),
        d_packet_length(0),
        d_frame_length(0),
        d_packet_num(0),
        d_on_packet(false),
        d_len_tag_key(pmt::string_to_symbol("packet_length")),
        d_frame_len_tag_key(pmt::string_to_symbol("frame_length")),
        d_num_tag_key(pmt::string_to_symbol("packet_num"))
    {
      d_header_data = new unsigned char[header_formatter->header_len()/d_dim];
      // Check if header_length is a multiple of the constellation dimensionality
      if (d_header_length % d_dim != 0){
        throw std::invalid_argument("Header length must be a multiple of the constellation dimension.");
      }
      set_tag_propagation_policy(TPP_DONT);
      set_min_noutput_items(d_header_length);
    }

    /*
     * Our virtual destructor.
     */
    mimo_ofdm_header_reader_cc_impl::~mimo_ofdm_header_reader_cc_impl()
    {
    }

    void
    mimo_ofdm_header_reader_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    void
    mimo_ofdm_header_reader_cc_impl::demod_header(const gr_complex *src, unsigned char *dest){
      for (unsigned int j = 0; j < d_header_length / d_dim; ++j) {
        dest[j] = d_constellation->decision_maker(&src[j*d_dim]);
      }
    }

    bool
    mimo_ofdm_header_reader_cc_impl::parse_header() {
      std::vector<tag_t> header_tags;
      // Check header validity.
      if (!d_header_formatter->header_parser(const_cast<const unsigned char *>(d_header_data), header_tags)) {
        // Header is corrupted.
        return false;
      } else {
        // Valid header.
        // Read header.
        for (unsigned int c = 0; c < header_tags.size(); c++) {
          // Try the key in different keyholes.
          if (pmt::equal(header_tags[c].key, d_len_tag_key)){
            d_packet_length = pmt::to_long(header_tags[c].value);
          } else if (pmt::equal(header_tags[c].key, d_frame_len_tag_key)){
            d_frame_length = pmt::to_long(header_tags[c].value);
          } else if (pmt::equal(header_tags[c].key, d_num_tag_key)){
            d_packet_num = pmt::to_long(header_tags[c].value);
          } else{
            GR_LOG_INFO(d_logger, format("Unknown header tag %s")%pmt::symbol_to_string(header_tags[c].key));
          }
        }
        d_header_tags = header_tags;
        return true;
      }
    }

    void
    mimo_ofdm_header_reader_cc_impl::add_tags(uint32_t offset) {
      for (unsigned int c = 0; c < d_header_tags.size(); c++) {
        add_item_tag(0,
                     nitems_written(0)+offset,
                     d_header_tags[c].key,
                     d_header_tags[c].value);
      }
    }

    int
    mimo_ofdm_header_reader_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      gr_complex const *in = (const gr_complex*)input_items[0];
      unsigned char *out = (unsigned char*)output_items[0];
      uint32_t nconsumed = 0;
      uint32_t nwritten = 0;

      // Check if we are currently on a packet.
      if (d_on_packet){
        // We are at the beginning of the payload of the current packet, which is completely in this buffer.
        // Check if this packet is interrupted by a new packet.
        std::vector <gr::tag_t> interrupt_tags;
        get_tags_in_window(interrupt_tags, 0, 0, d_packet_length, d_key);
        if(interrupt_tags.size() > 0){
          // We get interrupted. Dump current packet to the beginning of the next packet.
          nconsumed = interrupt_tags[0].offset - nitems_read(0);
        } else {
          // We don't get interrupted.
          // Copy payload to output and add tags.
          memcpy(out, in, sizeof(gr_complex) * d_packet_length);
          add_tags(0);
          nconsumed = d_packet_length;
          nwritten = d_packet_length;
        }
        // We finished this packet.
        d_on_packet = false;
      }
      // We are not on a packet and search for the next packet.
      std::vector <gr::tag_t> tags;
      get_tags_in_window(tags, 0, nconsumed, noutput_items, d_key);
      // Dump everything until the next packet arrives.
      if(tags.size() == 0){
        // There are no tags in this buffer. Dump samples.
        nconsumed = noutput_items;
      } else {
        // There are tags in this buffer. Dump samples before the first tag.
        nconsumed = tags[0].offset - nitems_read(0);
      }
      // Iterate over tags.
      uint32_t segment_length;
      for (unsigned int i = 0; i < tags.size(); ++i){
        if (i < tags.size() - 1) {
          // This is not the last tag.
          segment_length = tags[i + 1].offset - tags[i].offset;
          // Check if the next tag arrives within this header.
          if (d_header_length > segment_length){
            // Dump this samples.
            nconsumed += segment_length;
            continue;
          }
          // Demodulate header.
          demod_header(&in[nconsumed], d_header_data);
          // Consume header.
          nconsumed += d_header_length;
          // Parse header.
          if(parse_header()){
            // This header is valid.
            // Check if this packet is interrupted by a new packet.
            if (d_packet_length > segment_length-d_header_length){
              // We get interrupted. Dump current packet to the beginning of the next packet.
              nconsumed += segment_length-d_header_length;
            } else{
              // We don't get interrupted.
              // Copy payload to output and add tags.
              memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex) * d_packet_length);
              add_tags(nwritten);
              nconsumed += segment_length-d_header_length;
              nwritten += d_packet_length;
            }
          } else {
            // This header is invalid.
            //GR_LOG_INFO(d_logger, format("Invalid header at %d.") %(nitems_read(0)+nconsumed-d_header_length));
            // Dump the segment.
            nconsumed += segment_length-d_header_length;
          }
        } else {
          // This is the last tag.
          segment_length = noutput_items - nconsumed;
          // Check if there are header_length samples left in the input buffer.
          if (d_header_length > segment_length) {
            // Process header in the next work() call where the whole packet is available.
            set_min_noutput_items(d_header_length);
            break;
          }
          // Demodulate header.
          demod_header(&in[nconsumed], d_header_data);
          // Consume header.
          nconsumed += d_header_length;
          // Parse header.
          if(parse_header()){
            // This header is valid.
            // Check if there are packet_length samples left in the input buffer.
            if (d_packet_length > segment_length-d_header_length){
              // Process packet in the next work() call where the whole packet is available.
              d_on_packet = true;
              set_min_noutput_items(d_packet_length);
            } else{
              // Copy payload to output and add tags.
              memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex) * d_packet_length);
              add_tags(nwritten);
              nconsumed += segment_length-d_header_length;
              nwritten += d_packet_length;
              // Reset min_noutput_items to default for this block.
              set_min_noutput_items(d_header_length);
            }
          } else{
            // This header is invalid.
            //GR_LOG_INFO(d_logger, format("Invalid header at %d.") %(nitems_read(0)+nconsumed-d_header_length));
            // Dump the segment.
            nconsumed += segment_length-d_header_length;
          }
        }
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (nconsumed);

      // Tell runtime system how many output items we produced.
      return nwritten;
    }

  } /* namespace digital */
} /* namespace gr */

