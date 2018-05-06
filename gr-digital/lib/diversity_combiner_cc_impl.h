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
/*! \brief Diversity combining block with selection combining and maximum-ratio combining mode.
 *
 * \param num_inputs Number of inputs ports.
 * \param vlen Vector length of the input and output items.
 * \param combining_technique Combining technique selection combining (0) or maximum-ratio combining (1).
 */
    class diversity_combiner_cc_impl : public diversity_combiner_cc
    {
     private:
      uint16_t d_num_inputs; /*!< Number of inputs ports. */
      uint16_t d_vlen; /*!< Vector length of the input and output items. */
      uint8_t d_combining_technique;
      /*!< Combining technique selection combining (0) or maximum-ratio combining (1). */
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
      void process_symbol(gr_vector_const_void_star input, gr_complex* out, uint16_t offset, uint16_t length);
     public:
      diversity_combiner_cc_impl(uint16_t num_inputs, uint16_t vlen, uint8_t combining_technique);
      ~diversity_combiner_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_IMPL_H */

