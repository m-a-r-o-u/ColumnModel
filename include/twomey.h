#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <random>
#include "analize_state.h"
#include "grid.h"
#include "member_iterator.h"
#include "state.h"
#include "superparticle.h"
#include "superparticle_source.h"
#include "thermodynamic.h"
#include "ns_table.h"
#include "twomey_utils.h"

template <typename G>
double place_vertically_random(G& gen, State& state, int index) {
    std::uniform_real_distribution<> dis(-1, 1);
    double z = state.grid.getlay(index) + state.grid.length / 2. * dis(gen);
    if (!(z > state.grid.getlvl(index) && z < state.grid.getlvl(index + 1))) {
        throw std::out_of_range("particle with height " + std::to_string(z) + 
                                " is placed in the wrong layer by the placer_vertically_random routine");
    }
    return z;
}

inline double place_vertically_center(State& state, int index) {
    double z = state.grid.getlay(index);
    return z;
}
inline void feedback_qc(std::vector<Superparticle>& nuc_par, Layer& lay) {
    double qc_sum =
        std::accumulate(member_iterator(nuc_par.begin(), &Superparticle::qc),
                        member_iterator(nuc_par.end(), &Superparticle::qc), 0);
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
        : N_multi(0),
          N_sp(N_sp),
          nprf_cmp(N_lay, 0),
          Stab(N_sp, 0.),
          gen(gen) {
        Stab = nstable(N_sp, N_multi);
    };

    void init(Logger& logger) {
        logger.setAttr("N_sp", N_sp);
        logger.setAttr("N_multi", N_multi);
    }

    void generateParticles(OIt sp_itr, State& state, double dt,
                           const std::vector<Superparticle>& sp) {
        std::vector<int> nucprf = count_nucleated_ccn(sp, state.grid);
        std::vector<double> Sprf = supersaturation_profile(state);
        std::vector<int> nprf = indexes(Stab, Sprf);
        std::transform(nprf.begin(), nprf.end(), nprf.begin(), [this](int x){return x * this->N_multi;});

        for (auto itr = nprf.begin(); itr != nprf.end(); ++itr) {
            auto index = std::distance(nprf.begin(), itr);
            auto nuc_par = n_particle(*itr, nucprf[index], state, index);
            feedback_qc(nuc_par, state.layers[index]);
            append_par(nuc_par, sp_itr);
        }
    }

   private:
    std::vector<Superparticle> n_particle(int n, int n_cmp, State& state,
                                          int index) {
        int n_nuc = std::floor((n - n_cmp) / N_multi);
        if(n_nuc > N_sp){
        throw std::out_of_range("more particles will nucleate then the maximal amount: " + std::to_string(n_nuc) + 
                                "compare with (max) N_sp" + std::to_string(N_sp));
        }
        std::vector<Superparticle> res;
        if (n_nuc > 0) {
            res.reserve(n_nuc);
            for (int i = 0; i < n_nuc; ++i) {
                double r_crit = 1.e-7;
                double r_dry = std::min(r_crit, 8.e-10 / Stab[Stab.size()-1]);
                double r_init = std::min(r_crit, 8.e-10 / Stab[i]);
                int N = N_multi;
                double qc = cloud_water(N, r_init, r_dry, 1.);
                double z = place_vertically_random(gen, state, index);
                if (r_init < r_dry) {
                    throw std::logic_error("the inital radius r_init: " + std::to_string(r_init) + 
                                           "is larger then r_dry: " + std::to_string(r_dry));
                }
                res.push_back({qc, z, r_dry, N});
            }
            return res;
        } else {
            return res;
        }
    }

    int N_multi;
    const int N_sp;
    std::vector<int> nprf_cmp;
    std::vector<double> Stab;
    G& gen;
};

template <typename OIt, typename G>
class NoParticleSource : public SuperParticleSource<OIt> {
   public:
    NoParticleSource(G& gen, int N_sp, int N_lay)
        : N_multi(1.e8 / double(N_sp)), N_sp(N_sp), gen(gen){};

    void init(Logger& logger) {
        logger.setAttr("N_sp", N_sp);
        logger.setAttr("N_multi", N_multi);
    }

    void generateParticles(OIt sp_itr, State& state, double dt,
                           const std::vector<Superparticle>& sp) {}

   private:
    const int N_multi;
    const int N_sp;
    G& gen;
};

template <typename OIt, typename G>
std::unique_ptr<Twomey<OIt, G>> mkTwomey(G& gen, const int& N_sp,
                                         const int& N_lay) {
    return std::make_unique<Twomey<OIt, G>>(gen, N_sp, N_lay);
}
