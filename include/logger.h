#pragma once
#include <iomanip>
#include <memory>
#include "state.h"
#include "superparticle.h"
#include "grid.h"
#include "thermodynamic.h"

class Logger {
   public:
    virtual void log(const State& state,
                     const std::vector<Superparticle>& superparticles,
                     const Grid& grid) const = 0;
    virtual void finalize() const = 0;
};

class StdoutLogger : public Logger {
   public:
    inline void log(const State& state,
                    const std::vector<Superparticle>& superparticles,
                    const Grid& grid) const override {
        std::cout << std::endl;
        std::cout << "State at " << state.t << "\n";
        std::cout << "     layer         z         E         p         T       "
                     " qv         S   (par) z        qc         r     r_dry\n";
        for (unsigned int i = 0; i < state.layers.size(); ++i) {
            std::cout << std::setprecision(3) << std::setw(10) << i;
            std::cout << std::setprecision(3) << std::setw(10)
                      << grid.getlay(i);
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].E;
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].p;
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].T;
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].qv;
            std::cout << std::setprecision(3) << std::setw(10)
                      << saturation(state.layers[i].T, state.layers[i].p, state.layers[i].qv);
            if (i < superparticles.size()) {
                std::cout << std::setprecision(3) << std::setw(10)
                          << superparticles[i].z;
                std::cout << std::setprecision(3) << std::setw(10)
                          << superparticles[i].qc;
                std::cout << std::setprecision(3) << std::setw(10)
                          << radius(superparticles[i].qc, superparticles[i].N, superparticles[i].r_dry);
                std::cout << std::setprecision(3) << std::setw(10)
                          << superparticles[i].r_dry;
            }
            std::cout << "\n";
        }
    }
    inline void finalize() const override {}
};

class NetcdfLogger : public Logger {};

inline std::unique_ptr<Logger> createLogger() {
    return std::make_unique<StdoutLogger>();
}
