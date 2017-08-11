#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "grid.h"
#include "member_iterator.h"
#include "state.h"
#include "superparticle.h"
#include "thermodynamic.h"

class LogisticGenerator {
   public:
    LogisticGenerator(int N_sp) : N_sp(N_sp){};
    double operator()() {
        ++count;
        return i_logistic_function();
    };

   private:
    double i_logistic_function() {
        double k = 1.3;
        double S0 = 0.1;
        return std::exp(std::log(count / (double(N_sp) - count)) / k +
                        std::log(S0)) /
               100.;
    };

    int N_sp;
    int count = 0;
};

inline std::vector<double> supersaturation_profile(const State& state) {
    std::vector<double> res(state.layers.size());
    std::transform(state.layers.begin(), state.layers.end(), res.begin(),
                   [](Layer l) { return saturation(l.T, l.p, l.qv); });
    return res;
}

template <typename IT>
int find_index(IT first, IT last, double value) {
    return std::distance(first, std::lower_bound(first, last, value));
}

class Twomey {
   public:
    Twomey(int N_sp, int N_multi) : g(N_sp), N_multi(N_multi) {
        S_lookup = createNSTable_logistic(N_sp);
    };
    inline void generateParticles(State& state, const Grid& grid,
                                  std::vector<Superparticle>& superparticles) {
        std::vector<double> S_state = supersaturation_profile(state);

        std::vector<int> N_profile;
        for (auto s : S_state) {
            N_profile.push_back(
                find_index(S_lookup.begin(), S_lookup.end(), s));
        }

        std::vector<double> lvls = grid.getlvls();

        for (auto itr = N_profile.begin(); itr != N_profile.end(); ++itr) {
            // superparticles.push_back()
            auto index = std::distance(N_profile.begin(), itr);
            auto nucleated_particles = first_n_particles(*itr, index, grid);

            superparticles.insert(superparticles.begin(),
                                  nucleated_particles.begin(),
                                  nucleated_particles.end());
        }
    }

   private:
    std::vector<double> createNSTable_logistic(const int N_sp) {
        std::vector<double> res(N_sp - 1);  // delete -1 and fix
        std::generate(res.begin(), res.end(), g);
        return res;
    }

    std::vector<Superparticle> first_n_particles(const int& n, const int& index,
                                                 const Grid& grid) {
        std::vector<Superparticle> res;
        res.reserve(n);
        for (int i = 0; i < n; ++i) {
            double S_i = S_lookup[i];
            double r_init = 8.e-10 / S_i;
            double z = grid.getlay(index);
            double r_dry = 0.;
            int N = N_multi;
            bool is_nucleated = true;
            double qc = cloud_water(N, r_init, r_dry, 1.);
            res.push_back({qc, z, r_dry, N, is_nucleated});
        }
        return res;
    }

    LogisticGenerator g;
    const int N_multi;
    std::vector<double> S_lookup;
};
