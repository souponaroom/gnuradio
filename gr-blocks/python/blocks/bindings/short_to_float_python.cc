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
/* BINDTOOL_HEADER_FILE(short_to_float.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(7514aa143e4ccafbdb5263071965ea36)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/blocks/short_to_float.h>
// pydoc.h is automatically generated in the build directory
#include <short_to_float_pydoc.h>

void bind_short_to_float(py::module& m)
{

    using short_to_float = ::gr::blocks::short_to_float;


    py::class_<short_to_float,
               gr::sync_block,
               gr::block,
               gr::basic_block,
               std::shared_ptr<short_to_float>>(m, "short_to_float", D(short_to_float))

        .def(py::init(&short_to_float::make),
             py::arg("vlen") = 1,
             py::arg("scale") = 1.,
             D(short_to_float, make))


        .def("scale", &short_to_float::scale, D(short_to_float, scale))


        .def("set_scale",
             &short_to_float::set_scale,
             py::arg("scale"),
             D(short_to_float, set_scale))

        ;
}
