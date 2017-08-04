#pragma once
#include <cmath>
#include <algorithm>
#include <vector>
#include "grid.h"
#include "superparticle.h"
#include "thermodynamic.h"

inline std::vector<double> count_superparticles(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);

        int index = std::distance(lvls.begin(), ubound) - 1;
        res[index] += 1;
    }
    return res;
}

inline std::vector<double> count_nucleated(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            res[index] += 1;
        }
    }
    return res;
}

inline std::vector<double> calculate_qc_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            res[index] += sp.qc;
        }
    }
    return res;
}

inline std::vector<double> calculate_effective_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> r2(grid.n_lay, 0);
    std::vector<double> r3(grid.n_lay, 0);
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            r2[index] += std::pow(radius(sp.qc, sp.N, sp.r_dry), 2);
            r3[index] += std::pow(radius(sp.qc, sp.N, sp.r_dry), 3);
        }
    }
    std::transform(r3.begin(), r3.end(), r2.begin(), res.begin(),
                   std::divides<void>());

    std::replace_if(res.begin(), res.end(),
                    [](const double& a) { return std::isnan(a); }, 0);
    return res;
}

inline std::vector<double> calculate_maximal_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            res[index] = std::max(radius(sp.qc, sp.N, sp.r_dry), res[index]);
        }
    }
    return res;
}

inline std::vector<double> calculate_minimal_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            res[index] = std::min(radius(sp.qc, sp.N, sp.r_dry), res[index]);
        }
    }
    return res;
}

inline std::vector<double> calculate_mean_radius_profile(
    const std::vector<Superparticle>& superparticles, const Grid& grid) {
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> count(grid.n_lay, 0);
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            count[index] += 1;
            res[index] += radius(sp.qc, sp.N, sp.r_dry);
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
    std::vector<double> lvls = grid.getlvls();
    std::vector<double> count(grid.n_lay, 0);
    std::vector<double> r2(grid.n_lay, 0);
    std::vector<double> mean(grid.n_lay, 0);
    std::vector<double> res(grid.n_lay, 0);

    for (auto sp : superparticles) {
        if (sp.is_nucleated) {
            auto ubound = std::upper_bound(lvls.begin(), lvls.end(), sp.z);
            int index = std::distance(lvls.begin(), ubound) - 1;
            count[index] += 1;
            r2[index] += std::pow(radius(sp.qc, sp.N, sp.r_dry), 2);
            mean[index] += radius(sp.qc, sp.N, sp.r_dry);
        }
    }
    std::transform(mean.begin(), mean.end(), count.begin(), mean.begin(),
                   std::divides<void>());

    std::replace_if(mean.begin(), mean.end(),
                    [](const double& a) { return std::isnan(a); }, 0);

    std::transform(r2.begin(), r2.end(), count.begin(),
                   res.begin(), std::divides<void>());

    std::replace_if(res.begin(), res.end(),
                    [](const double& a) { return std::isnan(a); }, 0);

    std::transform(res.begin(), res.end(),
                   res.begin(), [](double a){return std::sqrt(a);});

    std::transform(res.begin(), res.end(), mean.begin(),
                   res.begin(), std::minus<void>());

    return res;
}
