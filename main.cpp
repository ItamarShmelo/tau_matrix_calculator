#include "src/tau_matrix_monte_carlo.hpp"
#include "src/units.hpp"

int main(){
    Vector energy_groups_center = {
    7.363e-11,
    2.807e-10,
    7.2562e-10,
    1.6022e-9,
    2.3592e-9,
    3.5376e-9,
    5.5118e-9,
    9.143e-9,
    1.6667e-8,
    2.3616e-8,
    3.4863e-8,
    3.8712e-8,
    4.3125e-8,
    4.8203e-8,
    5.407e-8,
    6.0876e-8,
    6.8804e-8,
    7.8082e-8,
    8.8989e-8,
    1.0187e-7,
    1.1717e-7,
    1.3542e-7,
    1.5733e-7,
    1.8377e-7,
    2.1587e-7,
    2.5509e-7,
    3.0332e-7};
    
    Vector energy_groups_boundries = {
    1.3806e-16,
    1.2713e-10,
    4.2038e-10,
    1.1407e-9,
    1.8957e-9,
    2.7936e-9,
    4.2159e-9,
    6.6415e-9,
    1.2537e-8,
    1.9273e-8,
    3.0367e-8,
    3.6654e-8,
    4.0762e-8,
    4.5478e-8,
    5.0916e-8,
    5.7209e-8,
    6.4524e-8,
    7.3062e-8,
    8.3074e-8,
    9.4868e-8,
    1.0883e-7,
    1.2545e-7,
    1.4533e-7,
    1.6924e-7,
    1.9818e-7,
    2.3347e-7,
    2.7657e-7,
    1.6022e-5
    };

    auto tau_engine = tau_matrix_monte_carlo_engine(energy_groups_center, energy_groups_boundries, 200000, 0);

    // double const T = 2.0*units::me_c2 / units::k_boltz;
    double constexpr T = 10.0*units::kev_to_kelvin;
    Matrix m = tau_engine.generate_S_matrix(T);
}