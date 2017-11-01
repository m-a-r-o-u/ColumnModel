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

template <typename E>
class HallCollisionKernal {
   public:
    HallCollisionKernal(E efficiencies)
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
    E efficiencies;
};

struct SpMassTendencies {
    double dqc;
    double dN;
};

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
            : out(out), dt(dt), params(params) {
            pc = std::distance(first, last);
            assert(pc >= 2);
            csps.reserve(pc);
            for (auto it = first; it != last; ++it) {
                double r = it->radius();
                csps.push_back({r, size_t(std::distance(first, it)), double(it->N),
                                  params.sedimentation.fall_speed(r), it->qc});
            }
            std::sort(csps.begin(), csps.end());
        }

        void calculate() {
            size_t isp;
            for (size_t i = 0; i < pc; ++i) {
                isp = csps[i].i;
                out[isp].dN = weights(i);
                out[isp].dqc = 4. / 3. * PI * RHO_H2O * csps[i].N * mass(i);
            }
        }

       private:
        double weights(size_t i) {
            auto r = csps[i].r;
            double internal_collisions = -params.collision_kernal(r, r, 0) *
                                         0.5 * csps[i].N * (csps[i].N - 1);
            double external_collisions = 0;
            for (auto j = i + 1; j < pc; ++j) {
                auto R = csps[j].r;
                external_collisions -=
                    params.collision_kernal(r, R,
                                            csps[i].fs - csps[j].fs) *
                    csps[i].N * csps[j].N;
            }
            return dt * (internal_collisions + external_collisions);
        }

        double mass(size_t i) {
            auto ri = csps[i].r;
            double from_smaller = 0;
            for (size_t j = 0; j < i; ++j) {
                auto rj = csps[j].r;
                from_smaller += params.collision_kernal(
                                    rj, ri, csps[i].fs - csps[j].fs) *
                                csps[j].N * rj * rj * rj;
            }
            double from_larger = 0;
            for (size_t j = i + 1; j < pc; ++j) {
                auto rj = csps[j].r;
                from_larger -= params.collision_kernal(
                                   ri, rj, csps[i].fs - csps[j].fs) *
                               csps[j].N * ri * ri * ri;
            }
            return dt * (from_smaller + from_larger);
        }

        struct CollideSp {
            double r;
            size_t i;
            double N;
            double fs;
            double qc;
            bool operator<(const CollideSp& other) const { return r < other.r; }
        };

        std::vector<CollideSp> csps;
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
    BoxCollisions<HallCollisionKernal<Efficiencies>> bc(sedi,
                                          HallCollisionKernal<Efficiencies>({}));
    return std::make_unique<
        BoxCollisionAdapter<BoxCollisions<HallCollisionKernal<Efficiencies>>>>(bc);
}

inline std::unique_ptr<Collisions> mkNCS() {
    return std::make_unique<NoCollisions>();
}
