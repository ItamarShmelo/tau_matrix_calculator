#include <limits>
#include <ctime>
#include <cmath>

#include "tau_matrix_monte_carlo.hpp"
#include "units.hpp"

static double constexpr signaling_NaN = std::numeric_limits<double>::signaling_NaN();

tau_matrix_monte_carlo_engine::tau_matrix_monte_carlo_engine(Vector const energy_groups_center_, 
                                                             Vector const energy_groups_boundries_, 
                                                             std::size_t const num_of_samples_) :
                            energy_groups_center(energy_groups_center_),
                            energy_groups_boundries(energy_groups_boundries_),
                            num_energy_groups(energy_groups_center.size()),
                            tau_temp(num_energy_groups, Vector(num_energy_groups, signaling_NaN)),
                            num_of_samples(num_of_samples_), 
                            T(signaling_NaN),
                            theta(signaling_NaN),
                            sum_1_bt(signaling_NaN),
                            Sb(signaling_NaN),
                            sample_uniform_01(
                                boost::random::mt19937(static_cast<unsigned int>(std::time(0))),
                                boost::random::uniform_01<>()
                            ) {}

tau_matrix_monte_carlo_engine::tau_matrix_monte_carlo_engine(Vector const energy_groups_center_, 
                                                             Vector const energy_groups_boundries_, 
                                                             std::size_t const num_of_samples_, 
                                                             std::size_t const seed) :
                            energy_groups_center(energy_groups_center_),
                            energy_groups_boundries(energy_groups_boundries_),
                            num_energy_groups(energy_groups_center.size()),
                            tau_temp(num_energy_groups, Vector(num_energy_groups, signaling_NaN)),
                            num_of_samples(num_of_samples_), 
                            T(std::numeric_limits<double>::signaling_NaN()), 
                            sample_uniform_01(
                                boost::random::mt19937(static_cast<unsigned int>(seed)),
                                boost::random::uniform_01<>()
                            ) {}

double tau_matrix_monte_carlo_engine::sample_gamma(){
    double const r0Sb = sample_uniform_01()*Sb;
    
    double const r1 = sample_uniform_01();

    if(r0Sb <= 1.0){
        double const r2 = sample_uniform_01();
        double const r3 = sample_uniform_01();

        return 1.0 - theta*std::log(r1*r2*r3);
    }

    if(r0Sb <= sum_1_bt){
        double const r2 = sample_uniform_01();

        return 1.0 - theta*std::log(r1*r2);
    }

    return 1.0 - theta*std::log(r1);
}

void tau_matrix_monte_carlo_engine::set_temperature(double const temperature){
    T = temperature;
    theta = units::k_boltz * T / units::me_c2;
    sum_1_bt = 1.0 + 1.0 / theta;
    Sb = 1.0 + 1.0 / theta + 0.5/(theta*theta); 
}
