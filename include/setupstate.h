#pragma once
#include <vector>
#include "constants.h"

double exponential_qv(double z, double qv0, double zc) {
    return qv0 * std::exp(-z / zc);
}

double linear_temperature(double z, double T0) { return T0 - LAPSE_RATE_A * z; }

double hydrostatic_pressure(double z, double p0) { return p0 - G * z; }
