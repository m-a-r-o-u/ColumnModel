#pragma once
#include "thermodynamic.h"

struct Superparticle {
    double qc;
    double z;
    double r_dry;
    int N;
    bool is_nucleated;
    double S_prime = 0;
    double w_prime = 0;
    double radius() const;
};

inline double Superparticle::radius() const { return ::radius(qc, N, r_dry, 1.); }
