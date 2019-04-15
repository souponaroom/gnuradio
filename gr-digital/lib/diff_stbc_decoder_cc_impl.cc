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

    diff_stbc_decoder_cc::sptr
    diff_stbc_decoder_cc::make(float phase_offset, uint32_t vlen,
                               const std::string &start_key)
    {
      return gnuradio::get_initial_sptr
        (new diff_stbc_decoder_cc_impl(phase_offset, vlen, start_key));
    }

    /*
     * The private constructor
     */
    diff_stbc_decoder_cc_impl::diff_stbc_decoder_cc_impl(float phase_offset,
                                                         uint32_t vlen,
                                                         const std::string &start_key)
      : gr::block("diff_stbc_decoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)*vlen),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_vlen(vlen),
        d_key(pmt::string_to_symbol(start_key)),
        d_basis_vecs(std::vector<gr_complex>(2, std::polar((float)M_SQRT1_2, phase_offset)))
    {
      // Init predecessor with dummy sequence.
      d_predecessor = new gr_complex[2*vlen];
      for (unsigned int i = 0; i < vlen*2; ++i) {
        d_predecessor[i] = d_basis_vecs[0];
      }
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2*vlen);
      set_tag_propagation_policy(TPP_DONT);
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
      ninput_items_required[0] = noutput_items/d_vlen;
    }

    void
    diff_stbc_decoder_cc_impl::decode_sequences(const gr_complex *prev_seq,
                                                const gr_complex *seq,
                                                gr_complex *out,
                                                uint32_t length) {
      if (length > 0) {
        for (unsigned int i = 0; i < d_vlen; ++i) {
          // Calculate the dot product of the received input sequences with the new basis.
          gr_complex r_1 = seq[0*d_vlen + i] * std::conj(prev_seq[0*d_vlen + i]) + std::conj(seq[1*d_vlen + i]) * prev_seq[1*d_vlen + i];
          gr_complex r_2 = seq[0*d_vlen + i] * std::conj(prev_seq[1*d_vlen + i]) - std::conj(seq[1*d_vlen + i]) * prev_seq[0*d_vlen + i];
          // Calculate the decoded (but not normalized!) samples and write them to the output buffer.
          out[0*d_vlen + i] = d_basis_vecs[0] * r_1 - std::conj(d_basis_vecs[1]) * r_2;
          out[1*d_vlen + i] = d_basis_vecs[1] * r_1 + std::conj(d_basis_vecs[0]) * r_2;
        }
        // Recursively decode the remaining sequences of this block.
        decode_sequences(seq, &seq[2*d_vlen], &out[2*d_vlen], length - 2);
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
      get_tags_in_window(tags, 0, 0, noutput_items/d_vlen, d_key);
      // Check if the next tag is on an uneven position.
      for (unsigned int j = 0; j < tags.size(); ++j) {
        if (tags[j].offset%2 != 0){
          // This should be prevented by the system developer in most cases.
          GR_LOG_DEBUG(d_logger, format("Detected start tag on uneven position (tag[%d].offset = %d).\n "
                                        "This differential STBC scheme works on sequences of 2 samples. "
                                        "If you are not really sure what you are doing, "
                                        "you should only set 'start' tags on even sample positions.")
                                 %0 %tags[0].offset);
          break;
        }
      }
      if (tags.size() > 0) {
        if(tags[0].offset == nitems_read(0)){
          // Tag on 1st item.
          // Iterate over data blocks between taqs.
          for (unsigned int i = 0; i+1 < tags.size(); ++i) {
            // This is not the last tag in the buffer.
            /* Calculate input block length. The output block length is
             * one sequence shorter than the input block length
             * due to differential coding. */
            input_block_length = (tags[i+1].offset-(tags[i+1].offset%2)) - (tags[i].offset-(tags[i].offset%2));
            // Decode sequences of this block.
            if (input_block_length > 2) {
              decode_sequences(&in[nconsumed*d_vlen], &in[(nconsumed + 2)*d_vlen], &out[nproduced], input_block_length - 2);
              add_item_tag(0, nitems_written(0) + nproduced, d_key, pmt::from_long(0));
              //GR_LOG_DEBUG(d_logger, format("%d")%(nitems_written(0) + nproduced));
              nproduced += (input_block_length - 2)*d_vlen;
            }
            nconsumed += input_block_length;

          }
        } else{
          // There are items before the first tag.
          // Process samples before the first tag.
          input_block_length = (tags[0].offset-((tags[0].offset)%2)) - nitems_read(0);
          // Decode remaining sequences of the current block before the first tag.
          decode_sequences(d_predecessor, in, out, input_block_length);
          nconsumed = input_block_length;
          nproduced = input_block_length*d_vlen;
        }
      } else{
        // Process all samples, because there is no tag in this buffer.
        decode_sequences(d_predecessor, in, out, noutput_items);
        nconsumed = noutput_items/d_vlen;
        nproduced = noutput_items;
      }
      // Remember last sequence as predecessor of the next buffer.
      memcpy(d_predecessor, &in[noutput_items-2*d_vlen], sizeof(gr_complex)*2*d_vlen);

      // Tell runtime system how many input items we consumed.
      consume_each (nconsumed);

      // Tell runtime system how many output items we produced.
      return nproduced;
    }

  } /* namespace digital */
} /* namespace gr */

