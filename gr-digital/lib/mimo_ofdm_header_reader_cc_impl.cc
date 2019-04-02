/* -*- c++ -*- */
/* 
 * Copyright 2018,2019 Free Software Foundation, Inc.
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

    mimo_ofdm_header_reader_cc::sptr
    mimo_ofdm_header_reader_cc::make(constellation_sptr constellation,
                                     const gr::digital::packet_header_default::sptr &header_formatter,
                                     const std::string &start_key, const std::string &len_tag_key,
                                     const std::string &frame_len_tag_key, const std::string &num_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_header_reader_cc_impl(constellation, header_formatter, start_key,
                                             len_tag_key, frame_len_tag_key, num_tag_key));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_header_reader_cc_impl::mimo_ofdm_header_reader_cc_impl(
            constellation_sptr constellation,
            const gr::digital::packet_header_default::sptr &header_formatter,
            const std::string &start_key, const std::string &len_tag_key,
            const std::string &frame_len_tag_key, const std::string &num_tag_key)
      : gr::block("mimo_ofdm_header_reader_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_constellation(constellation),
        d_dim(constellation->dimensionality()),
        d_header_formatter(header_formatter),
        d_header_length(header_formatter->header_len()),
        d_start_key(pmt::string_to_symbol(start_key)),
        d_packet_length(0),
        d_frame_length(0),
        d_packet_num(0),
        d_on_packet(false),
        d_remaining_packet_len(0),
        d_len_tag_key(pmt::string_to_symbol(len_tag_key)),
        d_frame_len_tag_key(pmt::string_to_symbol(frame_len_tag_key)),
        d_num_tag_key(pmt::string_to_symbol(num_tag_key))
    {
      // Allocate space for demodulated header data.
      d_header_data = new unsigned char[header_formatter->header_len()/d_dim];
      // Check if header_length is a multiple of the constellation dimensionality.
      if (d_header_length % d_dim != 0){
        throw std::invalid_argument((format("Header length (%d) must be a multiple "
                                                    "of the constellation dimension (%d).")
                                     % d_header_length
                                     % d_dim).str());
      }
      /* Don't propagate any tags, because we add new tags (header info). */
      set_tag_propagation_policy(TPP_DONT);
      // The input buffer must be large enough to carry at least the samples of one header.
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

    uint32_t
    mimo_ofdm_header_reader_cc_impl::locate_next_tag(uint32_t offset,
                                                     uint32_t buffer_length) {
      // Find tags in buffer.
      std::vector <gr::tag_t> tags;
      get_tags_in_window(tags, 0, offset, buffer_length, d_start_key);

      if (tags.size() > 0){
        // Take next tag.
        return tags[0].offset - nitems_read(0);
      } else {
        // No further tags in this buffer.
        return buffer_length;
      }
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
            //GR_LOG_INFO(d_logger, format("Received unknown header component: %s.")%pmt::symbol_to_string(header_tags[c].key));
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
      gr_complex *out = (gr_complex*)output_items[0];
      uint32_t nconsumed = 0;
      uint32_t nwritten = 0;

      // State machine.
      if (d_on_packet){
        // We are still on the packet.
        // Check if this is the beginning of the packet.
        if (d_remaining_packet_len == d_packet_length){
          // Set tags.
          add_tags(0);
        }
        // Copy remaining packet.
        if (d_remaining_packet_len < noutput_items){
          // Remaining packet is within this buffer.
          // Copy remaining packet.
          memcpy(out, in, sizeof(gr_complex) * d_remaining_packet_len);
          nconsumed = d_remaining_packet_len;
          nwritten = d_remaining_packet_len;
          d_remaining_packet_len = 0;
          d_on_packet = false;

        } else {
          // The remaining packet length exceeds the current buffer.
          // Copy whole buffer.
          memcpy(out, in, sizeof(gr_complex) * noutput_items);
          nconsumed = noutput_items;
          nwritten = noutput_items;
          d_remaining_packet_len -= noutput_items;
        }
      } else {
        // We are not on a packet.
        // Find next 'start' tag in buffer.
        uint32_t next_start_pos;
        next_start_pos = locate_next_tag(0, noutput_items);
        if (next_start_pos < noutput_items){
          // The next start tag is within this buffer.
          // Dump samples until this next tag.
          nconsumed = next_start_pos;
          if (noutput_items-next_start_pos < d_header_length){
            // The remaining buffer size is smaller than the header length.
            // Process header in the next work() call where the whole header is available.
          } else {
            // The whole header is in this buffer.
            // Demodulate header.
            demod_header(&in[nconsumed], d_header_data);
            // Consume header.
            nconsumed += d_header_length;
            // Parse header.
            if (parse_header()) {
              // This header is valid.
              d_on_packet = true;
              d_remaining_packet_len = d_packet_length;
            } else {
              // This header is invalid.
              // Dump the segment.
              d_on_packet = false;
            }
          }
        } else {
          // There are no start tags in this buffer.
          // Drop samples.
          nconsumed = noutput_items;
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

