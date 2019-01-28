/* -*- c++ -*- */
/* 
 * Copyright 2018, 2019 Free Software Foundation, Inc.
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


#ifndef INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_H
#define INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_interpolator.h>

namespace gr {
  namespace digital {
    /*!
     * \brief Combines the input streams to one output stream to increase the SNR.
     * \ingroup digital
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

    class DIGITAL_API diversity_combiner_cc : virtual public gr::sync_interpolator
    {
     public:
      typedef boost::shared_ptr<diversity_combiner_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::diversity_combiner_cc.
       *
       * To avoid accidental use of raw pointers, digital::diversity_combiner_cc's
       * constructor is in a private implementation
       * class. digital::diversity_combiner_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint16_t num_inputs, uint16_t vlen, std::string combining_technique);
      
      virtual void set_combining_technique(std::string combining_technique) = 0;
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIVERSITY_COMBINER_CC_H */

