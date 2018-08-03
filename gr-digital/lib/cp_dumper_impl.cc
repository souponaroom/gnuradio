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
#include "cp_dumper_impl.h"

namespace gr {
  namespace digital {

    cp_dumper::sptr
    cp_dumper::make(int m, int n, int offset, int vlen)
    {
      return gnuradio::get_initial_sptr
        (new cp_dumper_impl(m, n, offset, vlen));
    }

    /*
     * The private constructor
     */
    cp_dumper_impl::cp_dumper_impl(int m, int n, int offset, int vlen)
      : gr::block("cp_dumper",
              gr::io_signature::make(1, 1, sizeof(gr_complex)*vlen),
              gr::io_signature::make(1, 1, sizeof(gr_complex)*vlen)),
        d_m(m), d_n(n), d_offset(offset), d_vlen(vlen)
    {
      set_output_multiple(n);
    }

    /*
     * Our virtual destructor.
     */
    cp_dumper_impl::~cp_dumper_impl()
    {
    }

    void
    cp_dumper_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = (noutput_items/d_m)*d_n;
    }

    int
    cp_dumper_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      for (int i = 0; i < noutput_items/d_m; ++i) {
        for (int j = 0; j < d_m; ++j) {
          for (int k = 0; k < d_vlen; ++k) {
            out[(i*d_m + j)*d_vlen + k] = in[(i*d_n + j + d_offset)*d_vlen + k];
          }

        }
      }
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each ((noutput_items/d_m)*d_n);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */

