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
        d_symbol_counter(0),
        d_header_length(header_formatter->header_len()),
        d_packet_length(0),
        d_frame_length(0),
        d_packet_num(0),
        d_on_frame(false),
        d_len_tag_key(pmt::string_to_symbol("packet_length")),
        d_frame_len_tag_key(pmt::string_to_symbol("frame_length")),
        d_num_tag_key(pmt::string_to_symbol("packet_num"))
    {
      d_header_data = new unsigned char[header_formatter->header_len()/d_dim];
      // Check if header_length is a multiple of the constellation dimensionality
      if (d_header_length % d_dim != 0){
        throw std::invalid_argument("Header length must be a multiple of the constellation dimension.");
      }
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
        GR_LOG_INFO(d_logger, format("Packet len = %d, Frame len = %d") %d_packet_length %d_frame_length);
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

      // Collect all tags of the input buffer in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);
      GR_LOG_INFO(d_logger, format("Found %d tags.") %tags.size());

      uint16_t segment_length; // Number of items in the current segment.

      if(tags.size() == 0) { // Input buffer includes no tags at all.
        if(d_on_frame){
          if (d_symbol_counter == 0){
            /* If we not yet wrote the first item of the packet, we still have
             * to add the tags from the header to the stream. */
            add_tags(0);
          }
          if(noutput_items < d_packet_length-d_symbol_counter){
            // The whole input buffer is part of the current frame. We copy it to the output.
            memcpy(out, in, sizeof(gr_complex)*noutput_items);
            d_symbol_counter += noutput_items;
            nwritten = noutput_items;
            nconsumed = noutput_items;
          } else {
            // The frame ends within this input buffer and there is no tag afterwards.
            // Copy the rest of the frame and reset frame state.
            memcpy(out, in, sizeof(gr_complex)*(d_packet_length-d_symbol_counter));
            nwritten = d_packet_length-d_symbol_counter;
            d_on_frame = false;
            // Dump the rest of the input buffer.
            nconsumed = noutput_items;
          }
        } else{
          // Dump complete input buffer, because we are not on a frame.
          nconsumed = noutput_items;
        }
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          // There are items in the input buffer, before the first tag arrives.
          segment_length = tags[0].offset - nitems_read(0);
          if(d_on_frame){
            if (d_symbol_counter == 0){
              /* If we not yet wrote the first item of the packet, we still have
               * to add the tags from the header to the stream. */
              add_tags(0);
            }
            if (segment_length >= d_packet_length-d_symbol_counter){
              // The frame ends within this input buffer, but there is no tag directly afterwards.
              // Copy the rest of the frame and reset frame state.
              memcpy(out, in, sizeof(gr_complex)*(d_packet_length-d_symbol_counter));
              nwritten = d_packet_length-d_symbol_counter;
              d_on_frame = false;
              // Dump the rest of the input buffer til the next tag position.
              nconsumed = segment_length;
            } else{
              // Copy data until the next tag.
              memcpy(out, in, sizeof(gr_complex)*segment_length);
              nwritten = segment_length;
              nconsumed = segment_length;
              d_on_frame = false;
            }
          } else{
            // Dump input til the next tag position, because we are not on a frame.
            nconsumed = segment_length;
          }
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i){
          // We are on a new frame.
          // Distinguish last tag from other tags.
          if (i < tags.size() - 1) {
            // This is not the last tag.
            // Calculate the segment length.
            segment_length = tags[i + 1].offset - tags[i].offset;
            // Check if whole header is in this buffer.
            if (segment_length < d_header_length){
              // The segment is shorter than the header length. Dump segment.
              nconsumed += segment_length;
              GR_LOG_INFO(d_logger, format("The distance between two start tags is shorter than the header length."));
              continue;
            }
            // Demodulate header.
            demod_header(&in[nconsumed], d_header_data);
            nconsumed += d_header_length;
            // Parse header.
            if(parse_header()){
              // The header is valid.
              GR_LOG_INFO(d_logger, format("Detected a valid packet at item %1%") % (nitems_read(0)+nconsumed));
              // Copy data until the end of the packet or the end of the segment is finished.
              memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex)*std::min(d_packet_length, segment_length-d_header_length));
              nwritten += std::min(d_packet_length, segment_length-d_header_length);
              nconsumed += segment_length-d_header_length;
              if (std::min(d_packet_length, segment_length-d_header_length) > 0){
                /* If we wrote at least one item to the output, we can add the tags
                 * from the header to the stream. */
                add_tags(nwritten);
              }
            } else{
              // The header is corrupted.
              GR_LOG_INFO(d_logger, format("Detected an invalid packet at item %1%") % (nitems_read(0)+nconsumed));
              // Dump the segment.
              nconsumed += segment_length-d_header_length;
            }
          } else {

            // This is the last tag.
            // Calculate the segment length.
            segment_length = noutput_items - (tags[i].offset - nitems_read(0));
            GR_LOG_INFO(d_logger, format("This is the last tag, segment length %d.")%segment_length);
            // Check if whole header is in this buffer.
            if (segment_length < d_header_length){
              /* The segment is shorter than the header length. But
               * it can continue with the next input buffer. Don't consume it
               * and process it with the next frame. */
              break;
            }
            // Demodulate header.
            demod_header(&in[nconsumed], d_header_data);
            nconsumed += d_header_length;
            // Parse header.
            if(parse_header()){
              // The header is valid.
              GR_LOG_INFO(d_logger, format("Detected a valid packet at item %1%") % (nitems_read(0)+nconsumed));
              // Check if the packet is finished within this buffer.
              if (d_packet_length <= segment_length-d_header_length){
                // The packet is finished within this buffer.
                // Copy packet and consume segment.
                memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex)*d_packet_length);
                nwritten += d_packet_length;
                nconsumed += segment_length-d_header_length;
                if (d_packet_length > 0){
                  /* If we wrote at least one item to the output, we can add the tags
                   * from the header to the stream. */
                  add_tags(nwritten);
                }
                d_on_frame = false;
              } else{
                // The packet continues after this buffer.
                // Copy segment (only part of the whole packet).
                memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex)*(segment_length-d_header_length));
                nwritten += (segment_length-d_header_length);
                nconsumed += segment_length-d_header_length;
                // Set packet counter to continue with this packet in the next general_work() call.
                d_on_frame = true;
                d_symbol_counter = segment_length-d_header_length;
                if (d_symbol_counter > 0){
                  /* If we wrote at least one item to the output, we can add the tags
                   * from the header to the stream. */
                  add_tags(nwritten);
                }
              }
            } else{
              // The header is corrupted.
              GR_LOG_INFO(d_logger, format("Detected an invalid packet at item %1%") % (nitems_read(0)+nconsumed));
              // Dump the segment.
              nconsumed += segment_length-d_header_length;
              d_on_frame = false;
            }
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

