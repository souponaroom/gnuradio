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
/* BINDTOOL_HEADER_FILE(prefs.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(f9d239804a24578ed1abfd735b31ca6b)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/prefs.h>
// pydoc.h is automatically generated in the build directory
#include <prefs_pydoc.h>

void bind_prefs(py::module& m)
{

    using prefs = ::gr::prefs;


    py::class_<prefs, std::shared_ptr<prefs>>(m, "prefs", D(prefs))

        .def(py::init<>(), D(prefs, prefs))


        .def("add_config_file",
             &prefs::add_config_file,
             py::arg("configfile"),
             D(prefs, add_config_file))


        .def("to_string", &prefs::to_string, D(prefs, to_string))


        .def("save", &prefs::save, D(prefs, save))


        .def(
            "has_section", &prefs::has_section, py::arg("section"), D(prefs, has_section))


        .def("has_option",
             &prefs::has_option,
             py::arg("section"),
             py::arg("option"),
             D(prefs, has_option))


        .def("get_string",
             &prefs::get_string,
             py::arg("section"),
             py::arg("option"),
             py::arg("default_val"),
             D(prefs, get_string))


        .def("set_string",
             &prefs::set_string,
             py::arg("section"),
             py::arg("option"),
             py::arg("val"),
             D(prefs, set_string))


        .def("get_bool",
             &prefs::get_bool,
             py::arg("section"),
             py::arg("option"),
             py::arg("default_val"),
             D(prefs, get_bool))


        .def("set_bool",
             &prefs::set_bool,
             py::arg("section"),
             py::arg("option"),
             py::arg("val"),
             D(prefs, set_bool))


        .def("get_long",
             &prefs::get_long,
             py::arg("section"),
             py::arg("option"),
             py::arg("default_val"),
             D(prefs, get_long))


        .def("set_long",
             &prefs::set_long,
             py::arg("section"),
             py::arg("option"),
             py::arg("val"),
             D(prefs, set_long))


        .def("get_double",
             &prefs::get_double,
             py::arg("section"),
             py::arg("option"),
             py::arg("default_val"),
             D(prefs, get_double))


        .def("set_double",
             &prefs::set_double,
             py::arg("section"),
             py::arg("option"),
             py::arg("val"),
             D(prefs, set_double))

        ;


    py::module m_thread = m.def_submodule("thread");
}
