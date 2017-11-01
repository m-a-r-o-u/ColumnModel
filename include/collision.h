#pragma once
#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include "constants.h"
#include "efficiencies.h"
#include "efficiencies.h"
#include "interpolate.h"
#include "member_iterator.h"
#include "sedimentation.h"
#include "superparticle.h"
#include "thermodynamic.h"

class HallCollisionKernal {
   public:
    HallCollisionKernal(Efficiencies efficiencies)
        : efficiencies(efficiencies) {}
    double operator()(double r, double R, double dfs) const {
        if (R <= 0.) {
            std::cout << "R in hall_collision_kernal is zero of smaller: " << R
                      << std::endl;
        }
        return PI * (R + r) * (R + r) * std::abs(dfs) *
               efficiencies.collision_efficiency(R * 1.e6, r / R);
    }

   private:
    Efficiencies efficiencies;
};

struct SpMassTendencies {
    double dqc;
    int dN;
};

inline std::ostream& operator<<(std::ostream& os, const SpMassTendencies& t) {
    os << t.dqc << " " << t.dN;
    return os;
}

class Collisions {
   public:
    virtual ~Collisions() {}
    virtual std::vector<SpMassTendencies> collide(
        const std::vector<Superparticle>& sps, const Grid& grid, double dt) = 0;
    virtual bool needs_sorted_superparticles() = 0;
};

static inline bool sp_zcmp(const Superparticle& sp, double z) {
    return sp.z < z;
}

template <typename CollisionKernal>
class BoxCollisions {
   public:
    BoxCollisions(const Sedimentation& sedimentation,
                  CollisionKernal collision_kernal)
        : sedimentation(sedimentation), collision_kernal(collision_kernal) {}
    template <typename SpIt, typename TIt>
    void collide(SpIt first, SpIt last, TIt out, double dt) {
        size_t pc = std::distance(first, last);
        if (pc < 2) {
            return;
        }
        Collider<SpIt, TIt> collider(first, last, out, dt, *this);
        collider.calculate();
    }

   private:
    const Sedimentation& sedimentation;
    CollisionKernal collision_kernal;

    template <typename SpIt, typename TIt>
    class Collider {
       public:
        Collider(SpIt first, SpIt last, TIt out, double dt,
                 const BoxCollisions& params)
            : sps(first), out(out), dt(dt), params(params) {
            pc = std::distance(first, last);
            assert(pc >= 2);
            ridx.reserve(pc);
            for (auto it = first; it != last; ++it) {
                ridx.emplace_back(it->radius(), std::distance(first, it));
            }
            std::sort(ridx.begin(), ridx.end());
            fall_speeds.reserve(pc);
            std::transform(ridx.begin(), ridx.end(), fall_speeds.begin(),
                           [&](const auto& ri) {
                               return params.sedimentation.fall_speed(ri.first);
                           });
        }

        void calculate() {
            size_t isp;
            for (size_t i = 0; i < pc - 1; ++i) {
                isp = ridx[i].second;
                out[isp].dN = weights(i);
                out[isp].dqc = 4. / 3. * PI * RHO_H2O *
                                   (sps[isp].N + out[isp].dN) * mass(i) /
                                   (1 - out[isp].dN / sps[isp].N) -
                               sps[isp].qc;
            }
            isp = ridx[pc - 1].second;
            out[isp].dqc = mass(pc - 1);
        }

       private:
        int weights(size_t i) {
            auto r = ridx[i].first;
            auto isp = ridx[i].second;
            double internal_collisions = -params.collision_kernal(r, r, 0) *
                                         0.5 * sps[isp].N * (sps[isp].N - 1);
            double external_collisions = 0;
            for (auto j = i + 1; j < pc; ++j) {
                auto R = ridx[j].first;
                auto iSP = ridx[j].second;
                external_collisions -=
                    params.collision_kernal(r, R,
                                            fall_speeds[i] - fall_speeds[j]) *
                    sps[isp].N * sps[iSP].N;
            }
            return std::floor(dt * (internal_collisions + external_collisions));
        }

        double mass(size_t i) {
            auto ri = ridx[i].first;
            // auto isp = ridx[i].second;
            double from_smaller = 0;
            for (size_t j = 0; j < i; ++j) {
                auto rj = ridx[j].first;
                auto jsp = ridx[j].second;
                from_smaller += params.collision_kernal(
                                    ri, rj, fall_speeds[i] - fall_speeds[j]) *
                                sps[jsp].N * rj * rj * rj;
            }
            double from_larger = 0;
            for (size_t j = i + 1; j < pc; ++j) {
                auto rj = ridx[j].first;
                auto jsp = ridx[j].second;
                from_larger -= params.collision_kernal(
                                   ri, rj, fall_speeds[i] - fall_speeds[j]) *
                               sps[jsp].N * ri * ri * ri;
            }
            return dt * (ri * ri * ri + from_smaller + from_larger);
        }

        std::vector<std::pair<double, size_t>> ridx;
        std::vector<double> fall_speeds;
        SpIt sps;
        TIt out;
        double dt;
        size_t pc;
        const BoxCollisions& params;
    };

    template <typename, typename>
    friend class Collider;
};

template <typename C>
class BoxCollisionAdapter : public Collisions {
   public:
    BoxCollisionAdapter(const C& boxcollider) : boxcollider(boxcollider) {}

    std::vector<SpMassTendencies> collide(const std::vector<Superparticle>& sps,
                                          const Grid& grid,
                                          double dt) override {
        std::vector<SpMassTendencies> tendencies(sps.size());
        const auto& lvls = grid.getlvls();
        if (lvls.empty()) {
            return tendencies;
        }
        auto it1 = std::lower_bound(sps.begin(), sps.end(), lvls[0], sp_zcmp);
        for (size_t i = 1; i < lvls.size(); ++i) {
            auto it2 = std::lower_bound(it1, sps.end(), lvls[i], sp_zcmp);
            auto tit1 = tendencies.begin() + std::distance(sps.begin(), it1);
            boxcollider.collide(it1, it2, tit1, dt);

            it1 = it2;
        }
        return tendencies;
    }
    bool needs_sorted_superparticles() override { return true; }

   private:
    C boxcollider;
};

class NoCollisions : public Collisions {
   public:
    NoCollisions() {}
    std::vector<SpMassTendencies> collide(const std::vector<Superparticle>& sps,
                                          const Grid& grid,
                                          double dt) override {
        return std::vector<SpMassTendencies>(sps.size());
    }
    bool needs_sorted_superparticles() override { return false; }
};

inline std::unique_ptr<Collisions> mkHCS(const Sedimentation& sedi) {
    BoxCollisions<HallCollisionKernal> bc(sedi,
                                          HallCollisionKernal(Efficiencies()));
    return std::make_unique<
        BoxCollisionAdapter<BoxCollisions<HallCollisionKernal>>>(bc);
}

inline std::unique_ptr<Collisions> mkNCS() {
    return std::make_unique<NoCollisions>();
}
