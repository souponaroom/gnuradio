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
from mimo import mimo_technique as mimo

class mimo_encoder_cc(gr.hier_block2):
    """
    Hierarchical MIMO encoder block of the following structure:
    -1 input port
    -encoder block with selected MIMO algorithm (default 'none' is no block at all)
    -M output ports
    """

    def __init__(self, mimo_technique, M=2, vlen=1):
        gr.hier_block2.__init__(self,
            "mimo_encoder_cc",
            gr.io_signature(1, 1, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(M, M, gr.sizeof_gr_complex))  # Output signature

        # Dictionary translating mimo algorithm keys into encoder blocks.
        mimo_algorithm = {mimo.ALAMOUTI.value : digital.alamouti_encoder_cc_make(vlen=vlen),
                          mimo.DIFF_ALAMOUTI.value : digital.diff_stbc_encoder_cc_make(block_len=vlen),
                          mimo.VBLAST_ZF.value : digital.vblast_encoder_cc_make(M),
                          mimo.VBLAST_MMSE.value: digital.vblast_encoder_cc_make(M),
                          mimo.RX_DIVERSITY_SC.value : digital.vblast_encoder_cc_make(M),
                          mimo.RX_DIVERSITY_MRC.value: digital.vblast_encoder_cc_make(M)}

        # Check for valid M.
        if M < 1:
            raise ValueError('MIMO block must have M >= 1 (M=%d) selected).' % M)
        # Check for valid MIMO algorithm.
        if mimo_technique.value not in mimo_algorithm:
            raise ValueError('MIMO algorithm %s unknown.' % (mimo_technique))
        # Check if M = 2 for Alamouti-like schemes.
        if M != 2 and mimo_technique == ('alamouti' or 'diff_stbc'):
            raise ValueError('For Alamouti-like schemes like %s, M must be 2.' % mimo_technique)

        # Connect everything.
        mimo_encoder = mimo_algorithm[mimo_technique.value]
        self.connect(self, mimo_encoder)
        for m in range(0, M):
            self.connect((mimo_encoder, m), (self, m))
