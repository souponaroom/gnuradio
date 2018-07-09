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
from gnuradio import blocks, digital

class mimo_decoder_cc(gr.hier_block2):
    """
    Hierarchical MIMO decoder block of the following structure:
    -N input ports
    -channel estimator which produces stream tags with estimated channel matrix
    -decoder block with selected MIMO algorithm (default 'none' is no block at all)
    -1 output port

    For mimo_technique='none' and N>1, the decoder block only takes the first input port as output port.
    An estimated CSI is always tagged to the output stream, if a training sequence exists
    (also for the mimo_technique='none' case).
    """
    def __init__(self, M, N, mimo_technique='none', training_sequence=[]):
        gr.hier_block2.__init__(self,
            "mimo_decoder_cc",
            gr.io_signature(N, N, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(1, 1, gr.sizeof_gr_complex)) # Output signature

        # Dictionary translating mimo algorithm keys into encoder blocks.
        mimo_algorithm = {'none': self,
                          'diversity_combining_SC' : digital.diversity_combiner_cc_make(N, 0, 'SC'),
                          'diversity_combining_MRC' : digital.diversity_combiner_cc_make(N, 0, 'MRC'),
                          'alamouti': digital.alamouti_decoder_cc_make(),
                          'diff_stbc': digital.diff_stbc_decoder_cc_make(),
                          'vblast': digital.vblast_decoder_cc_make(N, 'ZF')}

        if len(training_sequence) > 0:
            channel_est = digital.mimo_channel_estimator_cc_make(M, N, training_sequence)

        # Connect everything.
        if mimo_technique != 'none':
            if len(training_sequence) > 0:
                for n in range(0, N):
                    self.connect((self, n), channel_est, (mimo_algorithm, n))
            else:
                for n in range(0, N):
                    self.connect((self, n), (mimo_algorithm, n))
            self.connect(mimo_algorithm, self)
        else:
            if len(training_sequence) > 0:
                self.connect((self, 0), (channel_est, 0), (self, 0))
                for n in range(1, N):
                    self.connect((self, n), (channel_est, n), blocks.null_sink_make(gr.sizeof_gr_complex))
            else:
                self.connect(self, self)

