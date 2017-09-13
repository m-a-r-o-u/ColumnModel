#pragma once
#include <algorithm>
#include <memory>
#include <random>
#include <cmath>
#include <cstdlib>
#include "analize_state.h"
#include "grid.h"
#include "member_iterator.h"
#include "state.h"
#include "superparticle.h"
#include "superparticle_source.h"
#include "thermodynamic.h"

inline int twomey_ns(double S, int N_sp) {
    double k = 1.3;
    double S0 = 0.1;
    return N_sp / (1. + std::exp(-k * (std::log(S) - std::log(S0))));
}

inline double twomey_ns_inverse(int N, int N_sp) {
    double k = 1.3;
    double S0 = 0.1;
    return std::exp(std::log(double(N) / (N_sp - double(N))) / k +
                    std::log(S0));
}

class LogisticGenerator {
   public:
    LogisticGenerator(int N_sp) : N_sp(N_sp){};
    double operator()() {
        ++count;
        return i_logistic_function();
    };

   private:
    double i_logistic_function() { return twomey_ns_inverse(count, N_sp); };

    int N_sp;
    int count = 0;
};

template <typename IT>
int find_index(IT first, IT last, double value) {
    return std::distance(first, std::lower_bound(first, last, value));
}

template <typename G>
inline double place_vertically(G& gen, State& state, int index) {
    std::uniform_real_distribution<> dis(-1, 1);
    double z = state.grid.getlay(index) +  state.grid.length / 2. * dis(gen);
    if (!(z > state.grid.getlvl(index) && z < state.grid.getlvl(index + 1))) {
        std::exit(0);
    }
    return z;
}

inline std::vector<int> index_profile(const std::vector<double>& S_prf,
                                      const std::vector<double>& S_tab) {
    std::vector<int> N_prf;
    for (auto s : S_prf) {
        N_prf.push_back(find_index(S_tab.begin(), S_tab.end(), s));
    }
    return N_prf;
}

inline void feedback_qc(std::vector<Superparticle>& nuc_par, Layer& lay) {
    double qc_sum =
        std::accumulate(member_iterator(nuc_par.begin(),&Superparticle::qc),
                        member_iterator(nuc_par.end(),  &Superparticle::qc), 0);
    if (lay.qv > qc_sum) {
        lay.qv -= qc_sum;
    } else {
        lay.qv = 0.;
    }
}

template <typename OIt>
inline void append_par(const std::vector<Superparticle>& nuc_par, OIt sp_itr) {
    for (auto s : nuc_par) {
        *sp_itr++ = s;
    }
}

template <typename OIt, typename G>
class Twomey : public SuperParticleSource<OIt> {
   public:
    Twomey(G& gen, int N_sp, int N_lay)
        : g(N_sp),
          N_multi(1.e8 / double(N_sp)),
          N_sp(N_sp),
          N_pro_cmp(N_lay, 0),
          S_tab(N_sp, 0.),
          gen(gen)
    {
        S_tab = createNSTable_logistic(N_sp);
        std::transform(S_tab.begin(), S_tab.end(), S_tab.begin(),
                       [](double x) { return std::max(x, 8.* 1.e-3);});
    };

    void init(Logger& logger){
        logger.setAttr("N_sp", N_sp);
        logger.setAttr("N_multi", N_multi);
    }

    void generateParticles(OIt sp_itr, State& state, double dt,
                                  const std::vector<Superparticle>& sp) {
        std::vector<double> S_state = supersaturation_profile(state);
        std::vector<int> n_nuc = count_nucleated(sp, state.grid);
        std::vector<int> N_pro = index_profile(S_state, S_tab);

        for (auto itr = N_pro.begin(); itr != N_pro.end(); ++itr) {
            auto index = std::distance(N_pro.begin(), itr);
            auto nuc_par = n_particle(*itr, n_nuc[index], state, index);
            feedback_qc(nuc_par, state.layers[index]);
            append_par(nuc_par, sp_itr);
        }
    }

   private:
    std::vector<double> createNSTable_logistic(const int N_sp) {
        std::vector<double> res(N_sp, 0);
        std::generate(res.begin(), res.end(), g);
        return res;
    }

    std::vector<Superparticle> n_particle(int n, int n_cmp, State& state,
                                          int index) {
        int n_nuc = n - n_cmp;
        std::vector<Superparticle> res;
        if (n_nuc > 0) {
            res.reserve(n_nuc);
            for (int i = n_cmp; i < n_cmp + n_nuc; ++i) {
                double r_dry = 1.e-7;
                //double r_init = 8.e-10 / S_tab[i]; // is what comes from Grab. but leads to errors
                double r_init = r_dry;
                int N = N_multi;
                bool is_nucleated = true;
                if (r_init < r_dry){
                    std::cout << "\a" << std::endl;
                    std::exit(0);}
                double qc = cloud_water(N, r_init, r_dry, 1.);
                double z = place_vertically(gen, state, index);
                res.push_back({qc, z, r_dry, N, is_nucleated});
            }
            return res;
        } else {
            return res;
        }
    }

    LogisticGenerator g;
    const int N_multi;
    const int N_sp;
    std::vector<int> N_pro_cmp;
    std::vector<double> S_tab;
    G& gen;
};

template <typename OIt, typename G>
std::unique_ptr<Twomey<OIt, G>> mkTwomey(G& gen, const int& N_sp,
                                      const int& N_lay) {
    return std::make_unique<Twomey<OIt, G>>(gen, N_sp, N_lay);
}
