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

from gnuradio import gr, blocks
import digital_swig as digital
from mimo import mimo_technique as mimo

class mimo_decoder_cc(gr.hier_block2):
    """
    Hierarchical MIMO decoder block of the following structure:
    -N input ports
    -decoder block with selected MIMO algorithm (default 'none' is no block at all)
    -1 output port
    """
    def __init__(self, mimo_technique, N=2, vlen=1, csi_key="csi"):
        gr.hier_block2.__init__(self,
            "mimo_decoder_cc",
            gr.io_signature(N, N, gr.sizeof_gr_complex*vlen),  # Input signature
            gr.io_signature(1, 1, gr.sizeof_gr_complex))  # Output signature

        # Dictionary translating mimo algorithm keys into decoder blocks.
        mimo_algorithm = {mimo.RX_DIVERSITY_SC.value: digital.diversity_combiner_cc_make(num_inputs=N,
                                                                                  vlen=vlen,
                                                                                  combining_technique='SC'),
                          mimo.RX_DIVERSITY_MRC.value: digital.diversity_combiner_cc_make(num_inputs=N,
                                                                                    vlen=vlen,
                                                                                    combining_technique='MRC'),
                          mimo.ALAMOUTI.value: digital.alamouti_decoder_cc_make(vlen=vlen),
                          mimo.DIFF_ALAMOUTI.value: digital.diff_stbc_decoder_cc_make(vlen=vlen),
                          mimo.VBLAST_ZF.value: digital.vblast_decoder_cc_make(num_inputs=N,
                                                                    equalizer_type='ZF',
                                                                    vlen=vlen),
                          mimo.VBLAST_MMSE.value: digital.vblast_decoder_cc_make(num_inputs=N,
                                                                         equalizer_type='MMSE',
                                                                         vlen=vlen)
                          }

        # Check for valid N.
        if N < 1:
            raise ValueError('MIMO block must have N >= 1 (N=%d) selected).' % N)
        # Check for valid MIMO algorithm.
        if mimo_technique.value not in mimo_algorithm:
            raise ValueError('MIMO algorithm %s unknown.' % str(mimo_technique))
        # Check if N = 2 for Alamouti-like schemes.
        if N != 1 and mimo_technique.value is (mimo.ALAMOUTI.value or mimo.DIFF_ALAMOUTI.value):
            raise ValueError('For Alamouti-like schemes like %s, N must be 2.' % str(mimo_technique))

        # Connect everything.
        mimo_decoder = mimo_algorithm[mimo_technique.value]
        for i in range(0, N):
            self.connect((self, i), (mimo_decoder, i))

        self.connect(mimo_decoder, self)

