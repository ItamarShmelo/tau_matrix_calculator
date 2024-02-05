#include <pybind11/pybind11.h>
#include "../units.hpp"

namespace units {
    void bind_units(pybind11::module& m){
        namespace py = pybind11;

        m.attr("me") = py::float_(me);
        m.attr("c") = py::float_(c);
        m.attr("me_c2") = py::float_(me_c2);

        m.attr("k_boltz") = py::float_(k_boltz);

        m.attr("sigma_thomson") = py::float_(sigma_thomson);

        m.attr("ev") = py::float_(ev);
        m.attr("ev_to_kelvin") = py::float_(ev_to_kelvin);
        m.attr("kev_to_kelvin") = py::float_(kev_to_kelvin);
        
        m.attr("Navogadro") = py::float_(Navogadro);
    }
} 

PYBIND11_MODULE(_units, m){
    m.doc() = "units c++ module";
    units::bind_units(m);
}

