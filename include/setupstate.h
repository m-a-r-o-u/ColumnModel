#pragma once
#include <vector>
#include "constants.h"

std::vector<double> exponential_qv(std::vector<double> z, double qv0,
                                   double zc) {
    std::vector<double> qv;
    for (auto el : z) {
        qv.push_back(qv0 * std::exp(-el / zc));
    }
    return qv;
}

std::vector<double> linear_temperature(const std::vector<double>& z,
                                       double T0) {
    std::vector<double> T;
    for (auto el : z) {
        T.push_back(T0 - LAPSE_RATE_A * el);
    }
    return T;
}
std::vector<double> hydrostatic_pressure(const std::vector<double>& z,
                                         double p0) {
    std::vector<double> p;
    for (auto el : z) {
        p.push_back(p0 - G * el);
    }
    return p;
}

std::vector<double> get_z_levels(double dz, double n) {
    std::vector<double> z;
    for (int i = 0; i < n; ++i) {
        z.push_back(i * dz);
    }
    return z;
}
