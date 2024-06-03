/*
 * Copyright 2024 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(1)                                                        */
/* BINDTOOL_HEADER_FILE(diff_stbc_decoder_cc.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(86f757a6f7938d12f1b68c9fa5c9fee1)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/digital/diff_stbc_decoder_cc.h>
// pydoc.h is automatically generated in the build directory
#include <diff_stbc_decoder_cc_pydoc.h>

void bind_diff_stbc_decoder_cc(py::module& m)
{

    using diff_stbc_decoder_cc = ::gr::digital::diff_stbc_decoder_cc;


    py::class_<diff_stbc_decoder_cc,
               gr::block,
               gr::basic_block,
               std::shared_ptr<diff_stbc_decoder_cc>>(
        m, "diff_stbc_decoder_cc", D(diff_stbc_decoder_cc))

        .def(py::init(&diff_stbc_decoder_cc::make),
             py::arg("phase_offset") = 0.,
             py::arg("vlen") = 1,
             py::arg("start_key") = "start",
             D(diff_stbc_decoder_cc, make))


        ;
}
