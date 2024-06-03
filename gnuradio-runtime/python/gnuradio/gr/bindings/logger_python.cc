/*
 * Copyright 2020 Free Software Foundation, Inc.
 * Copyright 2021,2022 Marcus Müller
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
/* BINDTOOL_HEADER_FILE(logger.h)                                                  */
/* BINDTOOL_HEADER_FILE_HASH(b6745a64dce4f006d5bece226c75a01e)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/common.h>

namespace py = pybind11;

#include <gnuradio/logger.h>
#include <logger_pydoc.h>

void bind_logger(py::module& m)
{
    py::enum_<spdlog::level::level_enum>(m, "log_levels")
        // Values directly from spdlog/common.h
        .value("trace", spdlog::level::trace)
        .value("debug", spdlog::level::debug)
        .value("info", spdlog::level::info)
        .value("warn", spdlog::level::warn)
        .value("err", spdlog::level::err)
        .value("error", spdlog::level::err)
        .value("critical", spdlog::level::critical)
        .value("off", spdlog::level::off);

    using logger = gr::logger;

    py::class_<logger, std::shared_ptr<logger>>(m, "logger", D(logger))

        .def(py::init<std::string>(), py::arg("logger_name"), D(logger, logger))
        .def(py::init<gr::logger const&>(), py::arg("arg0"))

        .def("set_level",
             py::overload_cast<const std::string&>(&logger::set_level),
             py::arg("level"),
             D(logger, set_level))
        .def("set_level",
             py::overload_cast<const gr::log_level>(&logger::set_level),
             py::arg("level"),
             D(logger, set_level))
        .def("get_level",
             py::overload_cast<std::string&>(&logger::get_level, py::const_),
             py::arg("level"),
             D(logger, get_level))
        .def("get_string_level", &logger::get_string_level, D(logger, get_string_level))
        .def(
            "trace",
            [](logger& log, const std::string& msg) { log.trace("{:s}", msg); },
            py::arg("msg"),
            D(logger, trace))
        .def(
            "debug",
            [](logger& log, const std::string& msg) { log.debug("{:s}", msg); },
            py::arg("msg"),
            D(logger, debug))
        .def(
            "info",
            [](logger& log, const std::string& msg) { log.info("{:s}", msg); },
            py::arg("msg"),
            D(logger, info))
        .def(
            "notice",
            [](logger& log, const std::string& msg) { log.notice("{:s}", msg); },
            py::arg("msg"),
            D(logger, notice))
        .def(
            "warn",
            [](logger& log, const std::string& msg) { log.warn("{:s}", msg); },
            py::arg("msg"),
            D(logger, warn))
        .def(
            "error",
            [](logger& log, const std::string& msg) { log.error("{:s}", msg); },
            py::arg("msg"),
            D(logger, error))
        .def(
            "crit",
            [](logger& log, const std::string& msg) { log.crit("{:s}", msg); },
            py::arg("msg"),
            D(logger, crit))
        .def(
            "alert",
            [](logger& log, const std::string& msg) { log.alert("{:s}", msg); },
            py::arg("msg"),
            D(logger, alert))
        .def(
            "fatal",
            [](logger& log, const std::string& msg) { log.fatal("{:s}", msg); },
            py::arg("msg"),
            D(logger, fatal))
        .def(
            "emerg",
            [](logger& log, const std::string& msg) { log.emerg("{:s}", msg); },
            py::arg("msg"),
            D(logger, emerg))
        .def(
            "log",
            [](logger& log, spdlog::level::level_enum level, const std::string& msg) {
                log.log(level, "{:s}", msg);
            },
            py::arg("level"),
            py::arg("msg"),
            D(logger, log));

    using logging = gr::logging;

    py::class_<logging, std::unique_ptr<logging, py::nodelete>>(m, "logging")
        .def(py::init([]() {
                 return std::unique_ptr<logging, py::nodelete>(&logging::singleton());
             }),
             D(logging, singleton))
        .def("default_level", &logging::default_level, D(logging, default_level))
        .def("debug_level", &logging::debug_level, D(logging, debug_level))
        .def("set_default_level",
             &logging::set_default_level,
             D(logging, set_default_level))
        .def("set_debug_level", &logging::set_debug_level, D(logging, set_debug_level))
        .def("add_default_sink",
             &logging::add_default_sink,
             py::arg("sink"),
             D(logging, add_default_sink))
        .def("add_debug_sink",
             &logging::add_debug_sink,
             py::arg("sink"),
             D(logging, add_debug_sink))
        .def("add_default_console_sink",
             &logging::add_default_console_sink,
             D(logging, add_default_console_sink))
        .def("add_debug_console_sink",
             &logging::add_debug_console_sink,
             D(logging, add_debug_console_sink))
        .def_property_readonly_static(
            "default_pattern",
            [](py::object) { return logging::default_pattern; },
            D(logging, default_pattern));
}
