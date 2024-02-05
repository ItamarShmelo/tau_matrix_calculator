#ifndef UNITS
#define UNITS

// Everything is in cgs
namespace units {
static constexpr double me = 9.1093837e-28; // [gr] electron mass
static constexpr double c = 2.99792458000e10 ; // [cm/s] speed of light
static constexpr double me_c2 = me*c*c; // [gr*cm^2/s^2] electron rest energy

static constexpr double k_boltz = 1.3807e-16; // [gr * cm^2 / s^2 / K] boltzmann constant

static constexpr double sigma_thomson = 0.665e-25;

static constexpr double ev = 1.602176634e-12;
static constexpr double ev_to_kelvin = ev/k_boltz;
static constexpr double kev_to_kelvin = ev_to_kelvin*1e3;

static constexpr double Navogadro = 6.02214076e23;
}
#endif