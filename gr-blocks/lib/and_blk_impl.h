/* -*- c++ -*- */
/*
 * Copyright 2012,2018 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */


#ifndef AND_BLK_IMPL_H
#define AND_BLK_IMPL_H

#include <gnuradio/blocks/and_blk.h>

namespace gr {
namespace blocks {

template <class T>
class BLOCKS_API and_blk_impl : public and_blk<T>
{
    const size_t d_vlen;

public:
    and_blk_impl(size_t vlen);

    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items) override;
};

} /* namespace blocks */
} /* namespace gr */

#endif /* AND_BLK_IMPL_H */
