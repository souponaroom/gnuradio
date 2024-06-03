/*
 * Copyright 2020 Free Software Foundation, Inc.
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
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(throttle.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(9259d01362558f938bdbfec179b6e164)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/blocks/throttle.h>
// pydoc.h is automatically generated in the build directory
#include <throttle_pydoc.h>

void bind_throttle(py::module& m)
{

    using throttle = ::gr::blocks::throttle;


    py::class_<throttle,
               gr::sync_block,
               gr::block,
               gr::basic_block,
               std::shared_ptr<throttle>>(m, "throttle", D(throttle))

        .def(py::init(&throttle::make),
             py::arg("itemsize"),
             py::arg("samples_per_sec"),
             py::arg("ignore_tags") = true,
             py::arg("maximum_items_per_chunk") = 0,
             D(throttle, make))


        .def("set_sample_rate",
             &throttle::set_sample_rate,
             py::arg("rate"),
             D(throttle, set_sample_rate))


        .def("sample_rate", &throttle::sample_rate, D(throttle, sample_rate))

        ;
}
