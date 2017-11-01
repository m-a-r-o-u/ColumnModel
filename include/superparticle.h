#pragma once
#include <ostream>
#include "thermodynamic.h"

class Superparticle {
   public:
    Superparticle() = default;
    Superparticle(double qc, double z, double r_dry, int N, double v = 0,
                  double S_prime = 0, double w_prime = 0)
        : qc(qc),
          z(z),
          r_dry(r_dry),
          N(N),
          v(v),
          S_prime(S_prime),
          w_prime(w_prime) {
        update();
    }
    double qc;
    double z;
    double r_dry;
    int N;
    bool is_nucleated = true;
    double v = 0;
    double S_prime = 0;
    double w_prime = 0;
    inline double radius() const { return _radius; }
    void update() {
        _radius = ::radius(qc, N, r_dry, 1.);
        nucleation();
    }

   private:
    double _radius = 0;
    void nucleation() {
        if (qc < 0) {
            is_nucleated = false;
        }
        if (z < 0) {
            std::cout << "particle reached the ground, r:" << radius()
                      << std::endl;
            is_nucleated = false;
        }
        if (N < 0) {
            is_nucleated = false;
        }
    }
};

inline std::ostream& operator<<(std::ostream& os, const Superparticle& rhs) {
    os << rhs.qc << "qc " << rhs.z << "z " << rhs.r_dry << "r_dry " << rhs.N
       << "N " << rhs.is_nucleated << "is_nucleated ";
    return os;
}

struct IndexSuperparticle {
    unsigned int i;
    Superparticle sp;
};

/*
class Superparticle {
   public:
    Superparticle(double qc, double z, double r_dry, int N, bool is_nucleated,
                  double v = 0, double S_prime = 0, double w_prime = 0)
        : _qc(qc),
          z(z),
          _r_dry(r_dry),
          _N(N),
          is_nucleated(is_nucleated),
          v(v),
          S_prime(S_prime),
          w_prime(w_prime) {}

    double z;
    double v = 0;
    double S_prime = 0;
    double w_prime = 0;

   private:
    double _qc;
    double _r_dry;

   public:
    bool is_nucleated;

   private:
    int _N;
    bool _valid = false;

   public:
    inline double radius() const {
        if (!_valid) {
            _radius = ::radius(qc, N, r_dry, 1.);
        }
        _valid = true;
        return _radius;
    }
    inline double qc() const { return _qc; }
    inline void qc(double qc) const {
        _qc = qc;
        valid = false;
    }
    inline double r_dry() const { return _r_dry; }
    inline void r_dry(double r_dry) const {
        _r_dry = r_dry;
        valid = false;
    }
    inline int N() const { return _N; }
    inline void N(int N) const {
        _N = N;
        valid = false;
    }
};
*/
