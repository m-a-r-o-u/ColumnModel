#pragma once
#include <iostream>
#include <cmath>
#include "tendencies.h"
#include "constants.h"

double saturation(double T, double p, double qv);

double saturation_vapor(double T, double p);

double critical_saturation(double r_dry, double T);
/** \brief checks if particle nucleates in the current environment
 *
 * \param r_dry particle radius [mu m]
 * \param S Supersaturation [1]
 * \param T Temperature [K]
 * \returns true if the particle nucleates
 *
 * Done using ...
 */
bool will_nucleate(double r_dry, double S, double T);

Tendencies condensation(double qc, double N, double r_dry, double S, double T, double E,
                        double dt);

double radius(double qc, double N, double r_min = 0., double rho = RHO_AIR);

double cloud_water(double N, double r, double r_min = 0., double rho = RHO_AIR);

double saturation_pressure(const double T);

double condensation_solver(const double r_old, const double es, const double T,
                           const double S, const double E, const double dt);

double diffusional_growth(const double r_old, const double es, const double T,
                          const double S, const double E, const double dt);
