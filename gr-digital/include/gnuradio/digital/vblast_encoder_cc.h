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


#ifndef INCLUDED_DIGITAL_VBLAST_ENCODER_CC_H
#define INCLUDED_DIGITAL_VBLAST_ENCODER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_decimator.h>

namespace gr {
  namespace digital {

    /*!
     * \brief VBLAST encoder.
     * \ingroup digital
     * Encodes a serial input data stream into a vertical transmission vector
     * which is transmitted at once over multiple antennas.
     *
     */
    class DIGITAL_API vblast_encoder_cc : virtual public gr::sync_decimator
    {
     public:
      typedef boost::shared_ptr<vblast_encoder_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::vblast_encoder_cc.
       *
       * To avoid accidental use of raw pointers, digital::vblast_encoder_cc's
       * constructor is in a private implementation
       * class. digital::vblast_encoder_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint16_t num_outputs);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_VBLAST_ENCODER_CC_H */

