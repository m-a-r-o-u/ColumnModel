#pragma once
#include <iostream>
#include <cmath>
#include "tendencies.h"
#include "constants.h"
#include "layer_quantities.h"
#include "level_quantities.h"
#include "grid.h"

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

Tendencies condensation(double qc, double N, double r_dry, double S, double T,
                        double E, double dt);

double radius(double qc, double N, double r_min = 0., double rho = RHO_AIR);

double cloud_water(double N, double r, double r_min = 0., double rho = RHO_AIR);

double saturation_pressure(const double T);

double condensation_solver(const double r_old, const double es, const double T,
                           const double S, const double E, const double dt);

double diffusional_growth(const double r_old, const double es, const double T,
                          const double S, const double E, const double dt);

template <typename QIt, typename WIt>
void advect_first_order(QIt q_begin, QIt q_end, WIt w_begin, double gridlength, double dt) {
    auto q_lo = *q_begin;
    ++q_begin;
    auto q_cur_it = q_begin;
    ++q_begin;
    auto q_hi_it = q_begin;

    auto w_lo = *w_begin;
    ++w_begin;
    auto w_hi_it = w_begin;
    double scale = dt / gridlength;

    while (q_hi_it != q_end) {
        auto q_cur_old = *q_cur_it;
        auto w_hi = *w_hi_it;
        auto q_hi = *q_hi_it;

        if (w_lo < 0) {
            *q_cur_it += scale * w_lo * q_cur_old;
        } else {
            *q_cur_it += scale * w_lo * q_lo;
        }
        if (w_hi < 0) {
            *q_cur_it -= scale * w_hi * q_hi;
        } else {
            *q_cur_it -= scale * w_hi * q_cur_old;
        }
        if (std::abs(w_lo * dt / gridlength) > 1 or std::abs(w_hi * dt / gridlength) > 1) {
            throw std::range_error("The CFL crit. is broken");
        }

        q_lo = q_cur_old;
        ++q_cur_it;
        ++q_hi_it;
        w_lo = *w_hi_it;
        ++w_hi_it;
    }
}

double fall_speed(const double r);
