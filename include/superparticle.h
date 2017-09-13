#pragma once
#include <ostream>
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

inline std::ostream& operator<<(std::ostream& os, const Superparticle& rhs){
    os  << rhs.qc << "qc "
        << rhs.z << "z "
        << rhs.r_dry << "r_dry "
        << rhs.N << "N "
        << rhs.is_nucleated << "is_nucleated ";
    return os;
}
inline double Superparticle::radius() const { return ::radius(qc, N, r_dry, 1.); }
