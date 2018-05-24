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

#include <gnuradio/io_signature.h>
#include "diff_stbc_decoder_cc_impl.h"

#include <boost/format.hpp>
using namespace boost;

namespace gr {
  namespace digital {

    const std::string diff_stbc_decoder_cc_impl::s = "start";
    const pmt::pmt_t diff_stbc_decoder_cc_impl::d_key = pmt::string_to_symbol(s);

    diff_stbc_decoder_cc::sptr
    diff_stbc_decoder_cc::make(float phase_offset)
    {
      return gnuradio::get_initial_sptr
        (new diff_stbc_decoder_cc_impl(phase_offset));
    }

    /*
     * The private constructor
     */
    diff_stbc_decoder_cc_impl::diff_stbc_decoder_cc_impl(float phase_offset)
      : gr::block("diff_stbc_decoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_basis_vecs(std::vector<gr_complex>(2, std::polar((float)M_SQRT1_2, phase_offset)))
    {
      // Init predecessor with dummy sequence.
      d_predecessor[0] = d_predecessor[1] = d_basis_vecs[0];
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2);
    }

    /*
     * Our virtual destructor.
     */
    diff_stbc_decoder_cc_impl::~diff_stbc_decoder_cc_impl()
    {
    }

    void
    diff_stbc_decoder_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    void
    diff_stbc_decoder_cc_impl::decode_sequences(const gr_complex *prev_seq,
                                            const gr_complex *seq,
                                            gr_complex *out,
                                            uint32_t length) {
      if(length > 0) {
        // Calculate the dot product of the received input sequences.
        gr_complex r_1 = seq[0] * std::conj(prev_seq[0]) + std::conj(seq[1]) * prev_seq[1];
        gr_complex r_2 = seq[0] * std::conj(prev_seq[1]) - std::conj(seq[1]) * prev_seq[0];
        // Calculate the decoded (but not normalized) samples and write them to the output buffer.
        out[0] = d_basis_vecs[0] * r_1 - std::conj(d_basis_vecs[1]) * r_2;
        out[1] = d_basis_vecs[1] * r_1 + std::conj(d_basis_vecs[0]) * r_2;

        // Recursively decode the remaining sequences of this block.
        decode_sequences(seq, &seq[2], &out[2], length-2);
      }
    }

    int
    diff_stbc_decoder_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      uint32_t nconsumed = 0;
      uint32_t nproduced = 0;
      uint32_t input_block_length;

      // Collect all tags of the input buffer with key "start" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);
      // Check if the next tag is on an uneven position.
      for (unsigned int j = 0; j < tags.size(); ++j) {
        if (tags[j].offset%2 != 0){
          // This should be prevented by the system developer in most cases.
          GR_LOG_DEBUG(d_logger, format("Detected start tag on uneven position (tag[%d].offset = %d).\n "
                                        "This differenatial STBC scheme works on sequences of 2 samples. "
                                        "If you are not really sure what you are doing, "
                                        "you should only set 'start' tags on even sample positions.")
                                 %0 %tags[0].offset);
          break;
        }
      }
      //GR_LOG_DEBUG(d_logger, format("%d noutput_items, %d tags")%noutput_items %tags.size());

      if (tags.size() > 0) {
        //GR_LOG_DEBUG(d_logger, format("There are tags at pos %d.")%tags[0].offset);
        // Process samples before the first tag.
        input_block_length = (tags[0].offset-(tags[0].offset%2)) - nitems_read(0);
        // Decode remaining sequences of the current block before the first tag.
        decode_sequences(d_predecessor, in, out, input_block_length);
        nconsumed = input_block_length;
        nproduced = input_block_length;

        //GR_LOG_DEBUG(d_logger, format("before tags: nconsumed: %d, nproduced %d")%nconsumed %nproduced);

        // Iterate over blocks between taqs.
        for (unsigned int i = 0; i+1 < tags.size(); ++i) {
          // This is not the last tag in the buffer.
          /* Calculate input block length. The output block length is
           * one sequence shorter than the input block length
           * due to differential coding. */
          input_block_length = (tags[i+1].offset-(tags[i+1].offset%2)) - (tags[i].offset-(tags[i].offset%2));
          //GR_LOG_DEBUG(d_logger, format("in iteration %d, input block length %d")%i %input_block_length);
          // Decode sequences of this block.
          if (input_block_length > 2) {
            decode_sequences(&in[nconsumed], &in[nconsumed + 2], &out[nproduced], input_block_length - 2);
            nproduced += input_block_length - 2;
          }
          nconsumed += input_block_length;

        }
        //GR_LOG_DEBUG(d_logger, format("After tags: nconsumed: %d, nproduced %d")%nconsumed %nproduced);

        // Process samples after last tag in the buffer.
        input_block_length = noutput_items - nconsumed;
        // Decode remaining sequences of this buffer.
        if (input_block_length - 2 > 0) {
          // There is more than one sequence left.
          decode_sequences(&in[nconsumed], &in[nconsumed + 2], &out[nproduced], input_block_length - 2);
          nproduced += input_block_length - 2;
        }
        nconsumed += input_block_length;
        //GR_LOG_DEBUG(d_logger, format("end: nconsumed: %d, nproduced %d")%nconsumed %nproduced);
      } else{
        //GR_LOG_DEBUG(d_logger, format("There are no tags."));
        // Process all samples, because there is no tag in this buffer.
        decode_sequences(d_predecessor, in, out, noutput_items);
        nconsumed = noutput_items;
        nproduced = noutput_items;
        //GR_LOG_DEBUG(d_logger, format("no tags: nconsumed: %d, nproduced %d")%nconsumed %nproduced);
      }

      //GR_LOG_DEBUG(d_logger, format("nconsumed: %d, nproduced %d")%nconsumed %nproduced);

      // Remember last sequence as predecessor of the next buffer.
      d_predecessor[0] = in[noutput_items-2];
      d_predecessor[1] = in[noutput_items-1];

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (nconsumed);

      // Tell runtime system how many output items we produced.
      return nproduced;
    }

  } /* namespace digital */
} /* namespace gr */

