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


#ifndef INCLUDED_DIGITAL_VBLAST_DECODER_CC_H
#define INCLUDED_DIGITAL_VBLAST_DECODER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_interpolator.h>

namespace gr {
  namespace digital {

/*! \brief VBLAST decoder.
 * Decodes incoming MxM MIMO data by equalizing the received signals in order to extract the
 * M different transmission signals. For the equalization, you can choose between a
 * zero forcing (key='ZF') and a minimum mean squared error (key='MMSE') scheme. Both schemes require
 * channel state information and MMSE additionally needs SNR information.
 *
 * The CSI and SNR information is transported via stream tags with key='csi' or 'snr', respectively.
 * The CSI must be a 3-dimensional vector with the dimensions vlen, num_inputs, num_inputs. For vlen>1
 * there is provided a separate channel matrix for each vector element. Set vlen>1 if you use a
 * multicarrier system like OFDM (in this case it would be vlen=number of occupied sub-carriers).
 * To generate a proper CSI tag, use the pmt structure pmt_vector(pmt_vector(pmt_c32vector))).
 * Initially the CSI is set to 1.0 + 0j for all branches and are updated
 * with each incoming CSI. The SNR is initially set to 1.0e6; in this case is the MMSE equalizer equal
 * to the ZF equalizer. The CSI and SNR tags are processed separately and they can therefore arrive
 * at different positions of the stream and occur in different frequencies.
 *
 * For 1x1 and 2x2 MIMO schemes, the equalizer is calculated internally. For MxM schemes with M > 2,
 * the C++ template library Eigen is used for the linear algebra operations.
 */
    class DIGITAL_API vblast_decoder_cc : virtual public gr::sync_interpolator
    {
     public:
      typedef boost::shared_ptr<vblast_decoder_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::vblast_decoder_cc.
       *
       * To avoid accidental use of raw pointers, digital::vblast_decoder_cc's
       * constructor is in a private implementation
       * class. digital::vblast_decoder_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint16_t num_inputs, std::string equalizer_type, uint16_t vlen=1);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_VBLAST_DECODER_CC_H */

