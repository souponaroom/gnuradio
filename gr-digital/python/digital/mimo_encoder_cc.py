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
from gnuradio import blocks
import digital_swig as digital
import pmt

class mimo_encoder_cc(gr.hier_block2):
    """
    Hierarchical MIMO encoder block of the following structure:
    -1 input port
    -encoder block with selected MIMO algorithm (default 'none' is no block at all)
    -insertion of training_sequence as pilot, at each MIMO output respectively
    -M output ports

    For mimo_technique='none' and M>1, the encoder block copies the input data to each output stream.
    """

    def __init__(self, M=1, mimo_technique='none', length_tag_name='length', training_sequence=[[]]):
        gr.hier_block2.__init__(self,
            "mimo_encoder_cc",
            gr.io_signature(1, 1, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(M, M, gr.sizeof_gr_complex)) # Output signature

        # Dictionary translating mimo algorithm keys into encoder blocks.
        mimo_algorithm = {'none' : self,
                          'alamouti' : digital.alamouti_encoder_cc_make(),
                          'diff_stbc' : digital.diff_stbc_encoder_cc_make(),
                          'vblast' : digital.vblast_encoder_cc_make(M)}

        # Check if M = 2 for Alamouti-like schemes.
        if M != 2 and mimo_technique == ('alamouti' or 'diff_stbc'):
            log = gr.logger("MIMO_logger")
            log.set_level("INFO")
            log.debug("M is fixed to 2 for alamouti-like MIMO schemes. Setting M=2.")
            M = 2

        # Connect everything.
        if mimo_technique != 'none':
            mimo_encoder = mimo_algorithm[mimo_technique]
            self.connect(self, mimo_encoder)
        else:
            mimo_encoder = self

        if len(training_sequence) > 0:
            mux = []
            training_src = []
            for m in range(0, M):
                mux.append(blocks.tagged_stream_mux_make(itemsize=gr.sizeof_gr_complex,
                                                         lengthtagname=length_tag_name,
                                                         tag_preserve_head_pos=1))
                tag = [gr.tag_utils.python_to_tag((0,
                                            pmt.string_to_symbol(length_tag_name),
                                            pmt.from_long(len(training_sequence)),
                                            pmt.from_long(0)))]

                self.connect(blocks.vector_source_c_make(data=training_sequence[m],
                                                         repeat=True,
                                                         vlen=1,
                                                         tags=tag),
                             (mux[m], 0))
                self.connect((mimo_encoder, m), (mux[m], 1), (self, m))
        else:
            for m in range(0, M):
                self.connect(mimo_encoder, (self, m))
