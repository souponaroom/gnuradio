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
from gnuradio import digital

class mimo_encoder_cc(gr.hier_block2):
    """
    docstring for block mimo_encoder_cc
    """

    def __init__(self, M=1, mimo_technique='none', training_sequence=[]):
        gr.hier_block2.__init__(self,
            "mimo_encoder_cc",
            gr.io_signature(1, 1, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(M, M, gr.sizeof_gr_complex)) # Output signature

        mimo_algorithm = {'none' : self,
                          'alamouti' : digital.alamouti_encoder_cc_make()}

        # Connect input to MIMO algorithm
        if mimo_technique != 'none':
            mimo_encoder = mimo_algorithm[mimo_technique]
            self.connect(self, mimo_encoder)
        else:
            mimo_encoder = self

        if len(training_sequence) > 0:
            for m in range(0, M):
                self.connect(mimo_encoder, insert_pilot, (self, m))
        else:
            for m in range(0, M)
                self.connect(mimo_encoder, (self, m))
