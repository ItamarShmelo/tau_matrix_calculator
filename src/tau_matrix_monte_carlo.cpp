#include <limits>
#include <ctime>
#include <cmath>

#include "tau_matrix_monte_carlo.hpp"
#include "units.hpp"

static double constexpr signaling_NaN = std::numeric_limits<double>::signaling_NaN();

tau_matrix_monte_carlo_engine::tau_matrix_monte_carlo_engine(Vector const energy_groups_center_, 
                                                             Vector const energy_groups_boundries_,
                                                             std::size_t const num_of_samples_,
                                                             bool const force_detailed_balance_) :
                            energy_groups_center(energy_groups_center_),
                            energy_groups_boundries(energy_groups_boundries_),
                            num_energy_groups(energy_groups_center.size()),
                            tau_temp(num_energy_groups, Vector(num_energy_groups, signaling_NaN)),
                            S_temp(num_energy_groups, Vector(num_energy_groups, signaling_NaN)),
                            num_of_samples(num_of_samples_), 
                            T(signaling_NaN),
                            theta(signaling_NaN),
                            sum_1_bt(signaling_NaN),
                            Sb(signaling_NaN),
                            sample_uniform_01(
                                boost::random::mt19937(static_cast<unsigned int>(std::time(0))),
                                boost::random::uniform_01<>()
                            ),
                            force_detailed_balance(force_detailed_balance_),

tau_matrix_monte_carlo_engine::tau_matrix_monte_carlo_engine(Vector const energy_groups_center_, 
                                                             Vector const energy_groups_boundries_, 
                                                             std::size_t const num_of_samples_, 
                                                             bool const force_detailed_balance_,
                                                             std::size_t const seed) :
                            energy_groups_center(energy_groups_center_),
                            energy_groups_boundries(energy_groups_boundries_),
                            num_energy_groups(energy_groups_center.size()),
                            tau_temp(num_energy_groups, Vector(num_energy_groups, signaling_NaN)),
                            S_temp(num_energy_groups, Vector(num_energy_groups, signaling_NaN)),
                            num_of_samples(num_of_samples_), 
                            T(std::numeric_limits<double>::signaling_NaN()), 
                            sample_uniform_01(
                                boost::random::mt19937(static_cast<unsigned int>(seed)),
                                boost::random::uniform_01<>()
                            ),
                            force_detailed_balance(force_detailed_balance_),

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

Matrix tau_matrix_monte_carlo_engine::generate_S_matrix(double const temperature, bool const log_grid){
    set_temperature(temperature);

    for(std::size_t i = 0; i < num_energy_groups; ++i){
        for(std::size_t j=0; j < num_energy_groups; ++j){
            S_temp[i][j] = 0.0;
        }
    }

    double sum_beta = 0.0;
    std::vector<double> const Omega_0(3., 0.);
    std::vector<double> Omega_0_tag(3., 0.), Omega_e(3., 0.), Omega_tag(3., 0.), Omega_p_tag(3., 0.);
    std::vector<double> weight(num_energy_groups, 0.);
    for(std::size_t sample_i=0; sample_i < num_of_samples; ++sample_i){
        // step 1: sample electron velocity from a weighted Maxwell Juttner distribution 
        double const gamma = sample_gamma();
        
        // weight of sample
        double const beta = std::sqrt(1.0 - 1.0 / (gamma*gamma));
        sum_beta += beta;

        // step 2: sample mu_e
        double const mu_e = 1.0 - 2.0*sample_uniform_01();
        Omega_e[0] = std::sqrt(1. - mu_e*mu_e);
        Omega_e[2] = mu_e;
        
        // step 3:
        double const D0 = gamma * (1.0 - beta*mu_e);
        
        // step 4: the direction of the photon before the scattering in the rest frame of the electron
        double const mu_0_tag = 1. / D0 * (1. - gamma/(1.+gamma)*(D0+1.)*beta*mu_e);
        double const sin_0_tag = std::sqrt(1. - mu_0_tag*mu_0_tag);
        Omega_0_tag[0] = -sin_0_tag;
        Omega_0_tag[2] = mu_0_tag;

        // step 5: sample the scattering angle of the photon
        double const mu_p_tag = 1.0 - 2.0*sample_uniform_01();
        double const sin_p_tag = std::sqrt(1.0 - mu_p_tag*mu_p_tag);
        double const psi_p_tag = sample_uniform_01()*2.*M_PI;
        Omega_p_tag[0] = sin_p_tag * std::cos(psi_p_tag);
        Omega_p_tag[1] = sin_p_tag * std::sin(psi_p_tag);
        Omega_p_tag[2] = mu_p_tag;

        // step 6 : rotate Omega_p by -theta_0 to get the angle in the electorn frame
        Omega_tag[0] = Omega_0_tag[2]*Omega_p_tag[0] + Omega_0_tag[0]*Omega_p_tag[2];
        Omega_tag[1] = Omega_p_tag[1];
        Omega_tag[2] = -Omega_0_tag[0]*Omega_p_tag[0] + Omega_0_tag[2]*Omega_p_tag[2];

        // step 7 
        double const D_tag = gamma*(1. + beta*(Omega_tag[0]*Omega_e[0] + Omega_tag[2]*Omega_e[2]));
        
        // step 8: sample the energy groups 
        double const interp = sample_uniform_01();
        for(std::size_t g0=0; g0<num_energy_groups; ++g0){
            // step 8a: sample energy
            double const E0 = energy_groups_boundries[g0] + interp*(energy_groups_boundries[g0+1]-energy_groups_boundries[g0]);
            // weight of energy sample
            double const w_E0 = E0*E0*std::exp(-E0/(units::k_boltz*T));
            weight[g0] += w_E0;
            
            // step 8b: calculate E
            double const E0_tag = D0*E0;
            double const A = 1. / (1. + (1. - mu_p_tag)*E0_tag / units::me_c2);
            double const E_tag = A*E0_tag;
            double const E = D_tag*E_tag;

            // step 8c: find the out energy group
            auto g_iterator = std::lower_bound(energy_groups_boundries.begin(), energy_groups_boundries.end(), E);
            auto g = std::distance(energy_groups_boundries.begin(), g_iterator)-1; // gives the index of the energy group

            g = std::max(0L, g);
            g = std::min(static_cast<long>(energy_groups_center.size())-1, g);

            // step 8d: calcualte the cross section contribution
            double const sigma = 0.75 * D0/gamma * A*A*(A + 1./A - sin_p_tag*sin_p_tag)*w_E0*beta;
            
            // step 9e: force the average energy change to be the monte carlo calculated energy change
            if(g0 == g){
                S_temp[g0][g] += sigma;
            } else {
                double const factor = std::min(1.0, (E-E0)/(energy_groups_center[g]-energy_groups_center[g0]));
                S_temp[g0][g] += factor*sigma;
                S_temp[g0][g0] += (1.0 - factor)*sigma;
            }
        }    
    }

    // total weight
    double const beta_avg = sum_beta / num_of_samples;
    
    // multiply by sigma_thomson and normalization factors
    for(std::size_t g0=0; g0 < num_energy_groups; ++g0){
        for(std::size_t g=0; g < num_energy_groups; ++g){
            double const weight_avg = weight[g0]/num_of_samples;
            S_temp[g0][g] *= units::sigma_thomson/(num_of_samples*beta_avg*weight_avg);
            if(log_grid) S_temp[g0][g] = std::log(S_temp[g0][g]);
        }
    }

    // Force detailed balance
    if(force_detailed_balance) detailed_balance();

    return S_temp;
}

void tau_matrix_monte_carlo_engine::detailed_balance(){
    double const k_bT = units::k_boltz*T;
    double constexpr thresh = units::sigma_thomson*std::numeric_limits<double>::epsilon()*1e3;
    for(std::size_t g0=0; g0<num_energy_groups; ++g0){
        double const w_g0 = energy_groups_boundries[g0+1]-energy_groups_boundries[g0];
        double const Eg0 = energy_groups_center[g0];
        double const Wien_g0 = Eg0*Eg0*std::exp(-Eg0/k_bT);
        for(std::size_t g=g0+1; g<num_energy_groups; ++g){
            if(S_temp[g0][g] < thresh && S_temp[g][g0] < thresh) continue;
            double const w_g = energy_groups_boundries[g+1]-energy_groups_boundries[g];
            double const Eg = energy_groups_center[g];
            double const Wien_g = Eg*Eg*std::exp(-Eg/k_bT);

            double const factor = w_g0*Wien_g0/(w_g*Wien_g);
            if(factor < 1.){
                S_temp[g][g0] = factor*S_temp[g0][g];
            } else {
                S_temp[g0][g] = factor*S_temp[g][g0];
            }
        }
    }
}

void tau_matrix_monte_carlo_engine::generate_S_log_tables(std::vector<double> const& tmp_grid){
    temperature_grid = tmp_grid;
    S_log_tables = std::vector<Matrix>(temperature_grid.size(), Matrix(num_energy_groups, Vector(num_energy_groups, 0.0)));

    for(std::size_t i=0; i < temperature_grid.size(); ++i){
        S_log_tables[i] = generate_S_matrix(temperature_grid[i], true);
    }
}
