#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../tau_matrix_monte_carlo.hpp"

namespace tau_matrix {
    void bind_tau_matrix_monte_carlo(pybind11::module& m){
        using namespace pybind11::literals;

        pybind11::class_<tau_matrix_monte_carlo_engine>(m, "tau_matrix_monte_carlo_engine")
        .def(pybind11::init<std::vector<double> const, 
                            std::vector<double> const, 
                            std::size_t const>(),
                            pybind11::kw_only(),
                            "energy_groups_center"_a,
                            "energy_groups_boundaries"_a,
                            "num_of_samples"_a)
        .def(pybind11::init<std::vector<double> const,
                            std::vector<double> const,
                            std::size_t const,
                            std::size_t const>(),
                            pybind11::kw_only(),
                            "energy_groups_center"_a,
                            "energy_groups_boundaries"_a,
                            "num_of_samples"_a,
                            "seed"_a)
        .def("sample_gamma",    &tau_matrix_monte_carlo_engine::sample_gamma)
        .def("set_temperature", &tau_matrix_monte_carlo_engine::set_temperature, pybind11::kw_only(),"temperature"_a)
        .def("generate_S_matrix", &tau_matrix_monte_carlo_engine::generate_S_matrix, pybind11::kw_only(), "temperature"_a)
        .def_readonly("theta",  &tau_matrix_monte_carlo_engine::theta)
        ;
    }
}


PYBIND11_MODULE(_tau_matrix_monte_carlo, m){
    m.doc() = "tau matrix monte carlo c++ module";

    tau_matrix::bind_tau_matrix_monte_carlo(m);
}