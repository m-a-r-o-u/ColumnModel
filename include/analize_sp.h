#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "grid.h"
#include "superparticle.h"
#include "thermodynamic.h"

inline void removeUnnucleated(std::vector<Superparticle>& superparticles) {
    auto fwd_it =
        std::remove_if(superparticles.begin(), superparticles.end(),
                       [](Superparticle s) { return !s.is_nucleated; });
    superparticles.erase(fwd_it, superparticles.end());
}

template <typename T>
inline std::vector<T> count_sp(
    const std::vector<Superparticle>& superparticles, 
    const Grid& grid,
    bool (*f)(const Superparticle& s), 
    T (*g)(const Superparticle& s))
{
    std::vector<T> res(grid.n_lay, 0);
    for (const auto& sp : superparticles) {
       if(f(sp)){
           int index =  grid.getlayindex(sp.z);
           res[index] += g(sp);
       }
    }
    return res;
}

inline std::vector<int> count_falling(const std::vector<Superparticle>& sps, const Grid& grid){
    return count_sp<int>(sps, grid, [](const Superparticle& s){return s.is_nucleated && (s.v<0);},
            [](const Superparticle& s){ return 1;});
}

inline std::vector<int> count_falling_ccn(const std::vector<Superparticle>& sps, const Grid& grid){
    return count_sp<int>(sps, grid, [](const Superparticle& s){return s.is_nucleated && (s.v<0);},
            [](const Superparticle& s){ return s.N;});
}

inline std::vector<int> count_nucleated(
    const std::vector<Superparticle>& sps, const Grid& grid) {
    return count_sp<int>(sps, grid, [](const Superparticle& s){return s.is_nucleated;},
            [](const Superparticle& s){return 1;});
}

inline std::vector<int> count_nucleated_ccn(
    const std::vector<Superparticle>& sps, const Grid& grid) {
    return count_sp<int>(sps, grid, [](const Superparticle& s){return s.is_nucleated;},
            [](const Superparticle& s){return s.N;});
}

inline std::vector<double> calculate_qc_profile(
    const std::vector<Superparticle>& sps, const Grid& grid) {
    return count_sp<double>(sps, grid, [](const Superparticle& s){return s.is_nucleated;},
            [](const Superparticle& s){return s.qc;});
}

inline std::vector<double> calculate_maximal_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> res(grid.n_lay, 0);
    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            int index =  grid.getlayindex(sp.z);
            res[index] = std::max(sp.radius(), res[index]);
        }
    }
    return res;
}


inline std::vector<double> calculate_minimal_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            int index =  grid.getlayindex(sp.z);
            res[index] = std::min(sp.radius(), res[index]);
        }
    }
    return res;
}


inline std::vector<double> calculate_effective_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> r2(grid.n_lay, 0);
    std::vector<double> r3(grid.n_lay, 0);
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            int index =  grid.getlayindex(sp.z);
            r2[index] += std::pow(sp.radius(), 2);
            r3[index] += std::pow(sp.radius(), 3);
        }
    }
    std::transform(r3.begin(), r3.end(), r2.begin(), res.begin(),
                   std::divides<void>());

    std::replace_if(res.begin(), res.end(),
                    [](const double& a) { return std::isnan(a); }, 0);
    return res;
}

inline std::vector<double> calculate_mean_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> count(grid.n_lay, 0);
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            int index =  grid.getlayindex(sp.z);
            count[index] += 1;
            res[index] += sp.radius();
        }
    }
    std::transform(res.begin(), res.end(), count.begin(), res.begin(),
                   std::divides<void>());

    std::replace_if(res.begin(), res.end(),
                    [](const double& a) { return std::isnan(a); }, 0);
    return res;
}

inline std::vector<double> calculate_stddev_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> count(grid.n_lay, 0);
    std::vector<double> r2(grid.n_lay, 0);
    std::vector<double> mean(grid.n_lay, 0);
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            int index =  grid.getlayindex(sp.z);
            count[index] += 1;
            r2[index] += std::pow(sp.radius(), 2);
            mean[index] += sp.radius();
        }
    }
    std::transform(mean.begin(), mean.end(), count.begin(), mean.begin(),
                   std::divides<void>());

    std::replace_if(mean.begin(), mean.end(),
                    [](const double& a) { return std::isnan(a); }, 0);

    std::transform(r2.begin(), r2.end(), count.begin(), res.begin(),
                   std::divides<void>());

    std::replace_if(res.begin(), res.end(),
                    [](const double& a) { return std::isnan(a); }, 0);

    std::transform(res.begin(), res.end(), res.begin(),
                   [](double a) { return std::sqrt(a); });

    std::transform(res.begin(), res.end(), mean.begin(), res.begin(),
                   std::minus<void>());

    return res;
}
