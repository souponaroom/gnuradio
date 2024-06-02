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
/* BINDTOOL_HEADER_FILE(vblast_encoder_cc.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(1fd6d8a0de8f7fdbc6f66efef81a7b66)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/digital/vblast_encoder_cc.h>
// pydoc.h is automatically generated in the build directory
#include <vblast_encoder_cc_pydoc.h>

void bind_vblast_encoder_cc(py::module& m)
{

    using vblast_encoder_cc = ::gr::digital::vblast_encoder_cc;


    py::class_<vblast_encoder_cc, gr::sync_decimator, std::shared_ptr<vblast_encoder_cc>>(
        m, "vblast_encoder_cc", D(vblast_encoder_cc))

        .def(py::init(&vblast_encoder_cc::make),
             py::arg("num_outputs"),
             D(vblast_encoder_cc, make))


        ;
}
