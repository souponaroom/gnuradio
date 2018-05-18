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

#ifndef INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_IMPL_H
#define INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_IMPL_H

#include <gnuradio/digital/diversity_combiner_cc.h>
#include <gnuradio/gr_complex.h>

namespace gr {
  namespace digital {
/*! \brief Combines the input streams to one output stream to increase the SNR.
 *
 * Diversity combining block. The following combining techniques are implemented right now:
 * - Selection Combining, key = 'SC'
 * - Maximum-Ratio Combining, key = 'MRC'
 *
 * The output stream is calculated out of a combination of the input streams with knowledge
 * of the channel state information (CSI) and with the algorithm of the selected combining technique.
 * The CSI is transported via stream tags with key='csi'. The combining parameters are initially
 * set to the selection of channel 0 for SC and an equal weighting of all channels for MRC
 * and are updated with each incoming CSI. The items between two tags are referred to as "symbols"
 * in this documentation. The length of a symbol can vary.
 *
 * \param num_inputs Number of inputs ports.
 * \param vlen Vector length of the input and output items.
 * \param combining_technique Combining technique. Selection combining ('SC') or
 * maximum-ratio combining ('MRC').
 */
    class diversity_combiner_cc_impl : public diversity_combiner_cc
    {
     private:
      uint16_t d_num_inputs; /*!< Number of inputs ports. */
      uint16_t d_vlen; /*!< Vector length of the input and output items. */
      std::string d_combining_technique;
      /*!< Combining technique selection combining ('SC') or maximum-ratio combining ('MRC'). */
      std::vector <gr::tag_t> tags; /*!< Vector that stores the tags in input buffer. */
      static const std::string s; /*!< String that matches the key of the CSI tags. */
      static const pmt::pmt_t d_key; /*!< PMT stores the key of the CSI tag. */
      std::vector<gr_complex> d_csi;
      /*!< Vector of length d_num_inputs which stores the current channel
       * state information (CSI). The vector is being updated which each
       * received tag of the key='csi'.
       */
      std::vector<float> d_csi_squared;
      /*!< Vector of length d_num_inputs which stores the current squared channel
       * state information (CSI). The vector is being updated which each
       * received tag of the key='csi'.
       */
      uint16_t d_best_path;
      /*!< Number of the input port which is selected as output for the current symbol. */
      std::vector<gr_complex> d_mrc_weighting;
      /*!< Vector of length d_num_inputs which stores the current normalized weighting vector. */

      void combine_inputs(gr_vector_const_void_star input, gr_complex* out, uint64_t offset, uint64_t length);
      /*! \brief Calculates the weighting vector out of CSI.
       *
       * @param input Vector of pointers to the input items, one entry per input stream.
       * @param out Pointer to the start of unwritten output items.
       * @param offset Number of already written input items.
       * @param length Number of items of the current symbol in the current buffer.
       */
      void process_symbol(gr_vector_const_void_star input, gr_complex* out, uint64_t offset, uint64_t length);
      /*! \brief Combines the input streams by the current weighting vector and
       * writes result to output buffer.
       *
       * @param input Vector of pointers to the input items, one entry per input stream.
       * @param out Pointer to the start of unwritten output items.
       * @param offset Number of already written input items.
       * @param length Number of items of the current symbol in the current buffer.
       */

    public:
      diversity_combiner_cc_impl(uint16_t num_inputs, uint16_t vlen, std::string combining_technique);
      ~diversity_combiner_cc_impl();

      void set_combining_technique(std::string combining_technique) {
        d_combining_technique = combining_technique;
      }

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_IMPL_H */

