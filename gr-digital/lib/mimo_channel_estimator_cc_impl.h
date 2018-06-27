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

#ifndef INCLUDED_DIGITAL_MIMO_CHANNEL_ESTIMATOR_CC_IMPL_H
#define INCLUDED_DIGITAL_MIMO_CHANNEL_ESTIMATOR_CC_IMPL_H

#include <gnuradio/digital/mimo_channel_estimator_cc.h>

namespace gr {
  namespace digital {

    class mimo_channel_estimator_cc_impl : public mimo_channel_estimator_cc
    {
     private:
      uint16_t d_num_inputs;
      std::vector<std::vector<gr_complex> > d_training_sequence;
      uint16_t d_training_length;
      std::vector <gr::tag_t> tags; /*!< Vector that stores the tags in input buffer. */
      static const pmt::pmt_t d_key; /*!< PMT stores the key of the CSI tag. */

      void copy_symbols(gr_vector_const_void_star &input_items,
                        gr_vector_void_star &output_items,
                        uint32_t symbol_length,
                        uint32_t reading_offset,
                        uint32_t writing_offset);

      void estimate_channel(gr_vector_const_void_star &input_items,
                            uint32_t reading_offset);

     public:
      mimo_channel_estimator_cc_impl(uint16_t num_inputs, std::vector<std::vector<gr_complex> > training_sequence);
      ~mimo_channel_estimator_cc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_MIMO_CHANNEL_ESTIMATOR_CC_IMPL_H */

