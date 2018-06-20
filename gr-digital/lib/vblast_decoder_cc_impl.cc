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

#include <string>
#include <boost/format.hpp>
#include <gnuradio/io_signature.h>
#include "vblast_decoder_cc_impl.h"
#include <eigen3/Eigen/Dense>

using namespace boost;

namespace gr {
  namespace digital {

    const std::string vblast_decoder_cc_impl::s = "csi";
    const pmt::pmt_t vblast_decoder_cc_impl::d_key = pmt::string_to_symbol(s);

    vblast_decoder_cc::sptr
    vblast_decoder_cc::make(uint16_t num_inputs, std::string equalizer_type)
    {
      return gnuradio::get_initial_sptr
        (new vblast_decoder_cc_impl(num_inputs, equalizer_type));
    }

    /*
     * The private constructor
     */
    vblast_decoder_cc_impl::vblast_decoder_cc_impl(uint16_t num_inputs, std::string equalizer_type)
      : gr::sync_interpolator("vblast_decoder_cc",
              gr::io_signature::make(num_inputs, num_inputs, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), num_inputs),
        d_num_inputs(num_inputs),
        d_equalizer_type(equalizer_type)
    {
      // Init CSI array and mimo_equalizer.
      d_csi = std::vector<std::vector<gr_complex> >(num_inputs, std::vector<gr_complex> (num_inputs, 1.0));
      d_mimo_equalizer = std::vector<std::vector<gr_complex> >(num_inputs, std::vector<gr_complex> (num_inputs, 1.0));
      d_snr = std::vector<float>(num_inputs, 1.0);
    }

    /*
     * Our virtual destructor.
     */
    vblast_decoder_cc_impl::~vblast_decoder_cc_impl()
    {
    }

    void
    vblast_decoder_cc_impl::update_mimo_equalizer() {
      if (d_equalizer_type.compare("ZF") == 0) {
        // Zero forcing equalizer.
        switch (d_num_inputs) {
          case 1: {
            // SISO case.
            d_mimo_equalizer[0][0] = (gr_complex) 1./d_csi[0][0];
            break;
          }
          case 2: {
            gr_complex c = d_csi[0][0] * d_csi[1][1] - d_csi[0][1] * d_csi[1][0];
            d_mimo_equalizer[0][0] = d_csi[1][1] / c;
            d_mimo_equalizer[0][1] = -d_csi[0][1] / c;
            d_mimo_equalizer[1][0] = -d_csi[1][0] / c;
            d_mimo_equalizer[1][1] = d_csi[0][0] / c;
            break;
          }
          default: {
            // Map CSI 2-dimensional std::vector to Eigen MatrixXcf.
            Eigen::MatrixXcf csi_matrix(d_num_inputs, d_num_inputs);
            for (int i = 0; i < d_num_inputs; ++i) {
              csi_matrix.row(i) = Eigen::VectorXcf::Map(&d_csi[i][0], d_csi[i].size());
            }
            // Calculate the inverse of the CSI matrix.
            Eigen::MatrixXcf csi_inverse = csi_matrix.inverse();

            // Map the inverse of the Eigen MatrixXcf to the 2-dim equalizer std::vector.
            for (int i = 0; i < d_num_inputs; ++i) {
              Eigen::VectorXcf::Map(&d_mimo_equalizer[i][0], csi_matrix.row(i).size()) = csi_inverse.row(i);
            }
          }
        }
      } else if (d_equalizer_type.compare("MMSE") == 0){
        // Minimum mean squared error equalizer.
        switch (d_num_inputs) {
          case 1: {
            // SISO case.
            d_mimo_equalizer[0][0] = std::conj(d_csi[0][0]) / (std::norm(d_csi[0][0])+(gr_complex) 1./d_snr[0]);
            break;
          }
          case 2: {
            gr_complex a = std::norm(d_csi[0][0]) + std::norm(d_csi[1][0]) + 1./d_snr[0];
            gr_complex b = std::conj(d_csi[0][0])*d_csi[0][1] + std::conj(d_csi[1][0])*d_csi[1][1];
            gr_complex c = std::conj(d_csi[0][1])*d_csi[0][0] + std::conj(d_csi[1][1])*d_csi[1][0];
            gr_complex d = std::norm(d_csi[0][1]) + std::norm(d_csi[1][1]) + 1./d_snr[1];
            gr_complex e = a*d-b*c;

            d_mimo_equalizer[0][0] = (std::conj(d_csi[0][0])*d - std::conj(d_csi[0][1])*b)/e;
            d_mimo_equalizer[0][1] = (std::conj(d_csi[1][0])*d - std::conj(d_csi[1][1])*b)/e;
            d_mimo_equalizer[1][0] = (std::conj(d_csi[0][1])*a - std::conj(d_csi[0][0])*c)/e;
            d_mimo_equalizer[1][1] = (std::conj(d_csi[1][1])*a - std::conj(d_csi[1][0])*c)/e;
            break;
          }
          default: {
            // Map CSI 2-dimensional std::vector to Eigen MatrixXcf.
            Eigen::MatrixXcf csi_matrix(d_num_inputs, d_num_inputs);
            for (int i = 0; i < d_num_inputs; ++i) {
              csi_matrix.row(i) = Eigen::VectorXcf::Map(&d_csi[i][0], d_csi[i].size());
            }
            // Calculate the pseudo-inverse fo the CSI matrix.
            Eigen::MatrixXcf snr_matrix = Eigen::MatrixXcf::Zero(d_num_inputs, d_num_inputs);
            for (int j = 0; j < d_num_inputs; ++j) {
              snr_matrix.coeffRef(j, j) = (gr_complex)(1.0 / d_snr[j]);
            }
            Eigen::MatrixXcf csi_pseudo_inverse = (csi_matrix.adjoint()*csi_matrix + snr_matrix).inverse() * csi_matrix.adjoint();

            // Map the inverse of the Eigen MatrixXcf to the 2-dim equalizer std::vector.
            for (int i = 0; i < d_num_inputs; ++i) {
              Eigen::VectorXcf::Map(&d_mimo_equalizer[i][0], csi_matrix.row(i).size()) = csi_pseudo_inverse.row(i);
            }
          }
        }
      } else{
        // The selected equalizer type is not existing in this implementation.
      }
    }

    void
    vblast_decoder_cc_impl::equalize_symbol(gr_vector_const_void_star input,
                                            gr_complex* out,
                                            uint32_t offset,
                                            uint32_t length) {
      std::fill(out, &out[length*d_num_inputs], 0.0);
      for (int n = 0; n < d_num_inputs; ++n) {
        gr_complex *in = &((gr_complex *) input[n])[offset/d_num_inputs];
        for (unsigned int i = 0; i < length; ++i) {
          for (int j = 0; j < d_num_inputs; ++j) {
            out[i*d_num_inputs + j] += d_mimo_equalizer[j][n] * in[i];
          }
        }
      }
    }

    int
    vblast_decoder_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      gr_complex *out = (gr_complex *) output_items[0];
      uint16_t nprocessed = 0; // Number of read and written items.

      // Collect all tags of the input buffer in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items);

      uint16_t symbol_length; // Number of items in the current symbol.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items;
        equalize_symbol(input_items, out, nprocessed, symbol_length);
        nprocessed += symbol_length;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
           * which belong to the previous symbol. */
          symbol_length = (tags[0].offset - nitems_read(0))*d_num_inputs;
          equalize_symbol(input_items, out, nprocessed, symbol_length);
          nprocessed += symbol_length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag.
            symbol_length = (tags[i + 1].offset - tags[i].offset)*d_num_inputs;
          } else {
            // This is the last tag.
            symbol_length = noutput_items - (tags[i].offset - nitems_read(0))*d_num_inputs;
          }
          // Check the key of the tag.
          if (pmt::symbol_to_string(tags[i].key).compare("csi") == 0) {
            // Get CSI from 'csi' tag.
            for (unsigned j = 0; j < pmt::length(tags[i].value); ++j) {
              d_csi[j] = pmt::c32vector_elements(pmt::vector_ref(tags[i].value, j));
            }
            // Calculate the weighting vector for the next symbol with the received CSI.
            update_mimo_equalizer();
          } else if (pmt::symbol_to_string(tags[i].key).compare("snr") == 0 && d_equalizer_type.compare("MMSE") == 0) {
            // 'snr' tag: Recalculate the weighting vector for the next symbol with the updated snr.
            d_snr = pmt::f32vector_elements(tags[i].value);
            update_mimo_equalizer();
          }
          // Process the symbol with the calculated weighting vector.
          equalize_symbol(input_items, &out[nprocessed], nprocessed, symbol_length);
          nprocessed += symbol_length;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */

