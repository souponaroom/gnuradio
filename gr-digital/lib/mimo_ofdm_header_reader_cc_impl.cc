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
        d_frame_length(0),
        d_on_frame(false)
    {
      d_header_data = new unsigned char[header_formatter->header_len()/d_dim];
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

      uint16_t segment_length; // Number of items in the current segment.

      if(tags.size() == 0) { // Input buffer includes no tags at all.
        if(d_on_frame){
          if(noutput_items < d_frame_length-d_symbol_counter){
            // The whole input buffer is part of the current frame. We copy it to the output.
            memcpy(out, in, sizeof(gr_complex)*noutput_items);
            d_symbol_counter += noutput_items;
            nwritten = noutput_items;
            nconsumed = noutput_items;
          } else {
            // The frame ends within this input buffer and there is no tag afterwards.
            // Copy the rest of the frame and reset frame state.
            memcpy(out, in, sizeof(gr_complex)*(d_frame_length-d_symbol_counter));
            nwritten = d_frame_length-d_symbol_counter;
            d_symbol_counter = 0;
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
            if (segment_length >= d_frame_length-d_symbol_counter){
              // The frame ends within this input buffer, but there is no tag directly afterwards.
              // Copy the rest of the frame and reset frame state.
              memcpy(out, in, sizeof(gr_complex)*(d_frame_length-d_symbol_counter));
              nwritten = d_frame_length-d_symbol_counter;
              d_symbol_counter = d_frame_length;
              // Dump the rest of the input buffer til the next tag position.
              nconsumed = segment_length;
            } else{
              // Copy data until the next tag.
              memcpy(out, in, sizeof(gr_complex)*segment_length);
              nwritten = segment_length;
              nconsumed = segment_length;
              d_symbol_counter += segment_length;
            }
          } else{
            // Dump input til the next tag position, because we are not on a frame.
            nconsumed = segment_length;
          }
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i){

          // Check if the last frame is unfinished.
          if (d_frame_length-d_symbol_counter > 0){
            /* The next tag comes before the current frame is finished.
               * We abort the current frame premature with a logger warning. */
            GR_LOG_WARN(d_logger, format("'start' tag before the end of the current frame. Aborting frame."));
          }
          // We are on a new frame. Reset symbol counter.
          d_symbol_counter = 0;
          // Calculate the segment length before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag.
            segment_length = tags[i + 1].offset - tags[i].offset;
          } else {
            // This is the last tag.
            segment_length = noutput_items - (tags[i].offset - nitems_read(0));
            d_on_frame = true;
          }
          // Read header TODO check if whole header symbol is in this buffer !!!
          for (int j = 0; j < d_header_formatter->header_len()/d_dim; ++j) { //TODO check if header_len/dim always int and if demod really to packed bits
            d_header_data[j] = d_constellation->decision_maker(&(in[nconsumed + j*d_dim]));
          }
          // Check header
          std::vector<tag_t> tags;
          if (!d_header_formatter->header_parser(const_cast<const unsigned char *>(d_header_data), tags)) {
            GR_LOG_INFO(d_logger, boost::format("Detected an invalid packet at item %1%") % (nitems_read(0)+nconsumed));
            d_frame_length = 0;
            d_on_frame = false;
          } else {
            // Valid header.
            d_frame_length = 0; // TODO read actual frame length.
            d_on_frame = false; // TODO change to true, if frame length read. !!!
            }
          nconsumed += d_header_formatter->header_len();
          // Process frame.
          if(segment_length >= d_frame_length-d_symbol_counter){
            // The frame ends within this input buffer, and there is no tag directly afterwards.
            // Copy the rest of the frame and reset frame state.
            memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex)*(d_frame_length-d_symbol_counter));
            nwritten = d_frame_length-d_symbol_counter;
            d_symbol_counter = 0;
            d_on_frame = false;
            // Dump the rest of the input buffer til the next tag position.
            nconsumed = segment_length;
          } else{
            // Copy data until the next tag.
            memcpy(&out[nwritten], &in[nconsumed], sizeof(gr_complex)*segment_length);
            nwritten = segment_length;
            nconsumed = segment_length;
            d_symbol_counter += segment_length;
          }
        }
      }
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */

