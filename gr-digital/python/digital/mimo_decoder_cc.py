#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2018 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr
import digital_swig as digital

class mimo_decoder_cc(gr.hier_block2):
    """
    Hierarchical MIMO decoder block of the following structure:
    -N input ports
    -decoder block with selected MIMO algorithm (default 'none' is no block at all)
    -1 output port
    """
    def __init__(self, N=2, mimo_technique='none', vlen=1):
        gr.hier_block2.__init__(self,
            "mimo_decoder_cc",
            gr.io_signature(N, N, gr.sizeof_gr_complex*vlen),  # Input signature
            gr.io_signature(1, 1, gr.sizeof_gr_complex))  # Output signature

        # Dictionary translating mimo algorithm keys into decoder blocks.
        mimo_algorithm = {'alamouti' : digital.alamouti_decoder_cc_make(),
                          'diff_stbc' : digital.diff_stbc_decoder_cc_make(),
                          'vblast' : digital.vblast_decoder_cc_make(num_inputs=N,
                                                                    equalizer_type='ZF',
                                                                    vlen=vlen)}

        # Check for valid N.
        if N < 1:
            raise ValueError('MIMO block must have N >= 1 (N=%d) selected).' % N)
        # Check for valid MIMO algorithm.
        if mimo_technique not in mimo_algorithm:
            raise ValueError('MIMO algorithm %s unknown.' % (mimo_technique))
        # Check if N = 2 for Alamouti-like schemes.
        if N != 2 and mimo_technique == ('alamouti' or 'diff_stbc'):
            raise ValueError('For Alamouti-like schemes like %s, N must be 2.' % mimo_technique)

        # Connect everything.
        mimo_decoder = mimo_algorithm[mimo_technique]
        for i in range(0, N):
            self.connect((self, i), (mimo_decoder, i))
        self.connect(mimo_decoder, self)

