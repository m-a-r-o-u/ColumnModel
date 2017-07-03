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
                     " qv         S     (lvl)         w   (par) z        qc\n";
        for (int i = 0; i < state.p.size(); ++i) {
            std::cout << std::setprecision(3) << std::setw(10) << i;
            std::cout << std::setprecision(3) << std::setw(10)
                      << grid.getlay(i);
            std::cout << std::setprecision(3) << std::setw(10) << state.E[i];
            std::cout << std::setprecision(3) << std::setw(10) << state.p[i];
            std::cout << std::setprecision(3) << std::setw(10) << state.T[i];
            std::cout << std::setprecision(3) << std::setw(10) << state.qv[i];
            std::cout << std::setprecision(3) << std::setw(10)
                      << saturation(state.T[i], state.p[i], state.qv[i]);
            std::cout << std::setprecision(3) << std::setw(10)
                      << grid.getlvl(i);
            std::cout << std::setprecision(3) << std::setw(10) << state.w[i];
            if (i < superparticles.size()) {
                std::cout << std::setprecision(3) << std::setw(10)
                          << superparticles[i].z;
                std::cout << std::setprecision(3) << std::setw(10)
                          << superparticles[i].qc;
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
