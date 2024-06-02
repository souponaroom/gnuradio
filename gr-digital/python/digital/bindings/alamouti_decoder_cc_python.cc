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
/* BINDTOOL_HEADER_FILE(alamouti_decoder_cc.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(573eb9db62b2f52c6098c795c6ac82eb)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/digital/alamouti_decoder_cc.h>
// pydoc.h is automatically generated in the build directory
#include <alamouti_decoder_cc_pydoc.h>

void bind_alamouti_decoder_cc(py::module& m)
{

    using alamouti_decoder_cc = ::gr::digital::alamouti_decoder_cc;


    py::class_<alamouti_decoder_cc,
               gr::sync_interpolator,
               std::shared_ptr<alamouti_decoder_cc>>(
        m, "alamouti_decoder_cc", D(alamouti_decoder_cc))

        .def(py::init(&alamouti_decoder_cc::make),
             py::arg("vlen") = 1,
             D(alamouti_decoder_cc, make))


        ;
}
