#pragma once
#include <iomanip>
#include <memory>
#include "state.h"
#include "superparticle.h"
#include "grid.h"
#include "thermodynamic.h"
#include "analize_distribution.h"

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
        std::vector<double> qc_sum = calculate_qc_profile(superparticles, grid);
        std::vector<double> r_mean = calculate_mean_radius_profile(superparticles, grid);
        std::vector<double> r_max = calculate_maximal_radius_profile(superparticles, grid);
        std::vector<double> sp_count = count_superparticles(superparticles, grid);
        std::vector<double> sp_count_nuc = count_nucleated(superparticles, grid);

        std::cout << std::endl;
        std::cout << "State at " << state.t << "\n";
        std::cout << "     layer         z         E         p         T       "
                     " qv         S        qc    r_mean     r_max     N_tot     N_nuc\n";
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
            std::cout << std::setprecision(3) << std::setw(10) << qc_sum[i];
            std::cout << std::setprecision(3) << std::setw(10) << r_mean[i];
            std::cout << std::setprecision(3) << std::setw(10) << r_max[i];
            std::cout << std::setprecision(3) << std::setw(10) << sp_count[i];
            std::cout << std::setprecision(3) << std::setw(10) << sp_count_nuc[i];
            std::cout << "\n";
        }
        std::cout << std::endl;
    }
    inline void finalize() const override {}
};

class NetcdfLogger : public Logger {};

inline std::unique_ptr<Logger> createLogger() {
    return std::make_unique<StdoutLogger>();
}
