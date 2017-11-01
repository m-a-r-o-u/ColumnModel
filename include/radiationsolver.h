#pragma once
#include <algorithm>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <iterator>
#include <string>
#include <vector>
#include "backgroundlevel.h"
#include "fpda_rrtm_lw_cld.h"
#include "fpda_rrtm_sw_cld.h"
#include "grid.h"
#include "readatm_utils.h"
#include "state.h"
#include "superparticle.h"

template <typename IIt, typename IIt2>
std::vector<double> concatonate(IIt a_first, IIt a_last, IIt2 b_first,
                                IIt2 b_last) {
    std::vector<double> res;
    res.reserve(a_last - a_first + b_last - b_first);
    res.insert(res.end(), a_first, a_last);
    res.insert(res.end(), b_first, b_last);
    return res;
}

static inline void calculate_cloudproperties(
    const std::vector<Superparticle>& superparticles, const Grid& grid,
    std::vector<double>& cliqwp, std::vector<double>& reliq) {
    std::vector<double> lvls = grid.getlvls();

    std::vector<double> qc_sum = calculate_qc_profile(superparticles, grid);
    std::vector<double> r_eff =
        calculate_effective_radius_profile(superparticles, grid);

    std::transform(qc_sum.begin(), qc_sum.end(), qc_sum.begin(),
                   std::bind(std::multiplies<void>(), std::placeholders::_1,
                             1.e3 * grid.length));

    std::transform(qc_sum.begin(), qc_sum.end(), cliqwp.begin(), cliqwp.begin(),
                   std::plus<void>());

    std::transform(
        r_eff.begin(), r_eff.end(), r_eff.begin(),
        std::bind(std::multiplies<void>(), std::placeholders::_1, 1.e6));

    std::transform(r_eff.begin(), r_eff.end(), reliq.begin(), reliq.begin(),
                   std::plus<void>());

    std::reverse(reliq.begin(), reliq.end());
    std::reverse(cliqwp.begin(), cliqwp.end());
    double re_min = 2.5;
    double re_max = 55.;
    std::replace_if(reliq.begin(), reliq.end(),
                    [re_min](double a) { return (a >= 0 && a < re_min); },
                    re_min);
    std::replace_if(reliq.begin(), reliq.end(),
                    [re_max](double a) { return (a > re_max); }, re_max);
}

struct RadiationSolver {
    RadiationSolver(std::string filename, bool sw, bool lw): sw(sw), lw(lw) {
        std::ifstream ifs(filename);
        std::vector<BackgroundLevelAfglus> bglvl;
        readin_atm<BackgroundLevelAfglus>(ifs, std::back_inserter(bglvl));

        z = get_quantity(bglvl, &BackgroundLevelAfglus::z);
        p = get_quantity(bglvl, &BackgroundLevelAfglus::p);
        T = get_quantity(bglvl, &BackgroundLevelAfglus::T);
        h2o = get_quantity(bglvl, &BackgroundLevelAfglus::h2o);
        o3 = get_quantity(bglvl, &BackgroundLevelAfglus::o3);
        o2 = get_quantity(bglvl, &BackgroundLevelAfglus::o2);
        no2 = get_quantity(bglvl, &BackgroundLevelAfglus::no2);
        co2 = get_quantity(bglvl, &BackgroundLevelAfglus::co2);
        air = get_quantity(bglvl, &BackgroundLevelAfglus::air);

        h2o = divide_elementwise(h2o, air);
        o3 = divide_elementwise(o3, air);
        o2 = divide_elementwise(o2, air);
        no2 = divide_elementwise(no2, air);
        co2 = divide_elementwise(co2, air);

        h2o = from_number_to_volume_ratio(h2o, M_MOL_H2O, RHO_H2O);
        o3 = from_number_to_volume_ratio(o3, M_MOL_O3, RHO_O3);
        o2 = from_number_to_volume_ratio(o2, M_MOL_O2, RHO_O2);
        no2 = from_number_to_volume_ratio(no2, M_MOL_NO2, RHO_NO2);
        co2 = from_number_to_volume_ratio(co2, M_MOL_CO2, RHO_CO2);

        T = pairwise_mean(T);
        h2o = pairwise_mean(h2o);
        o3 = pairwise_mean(o3);
        o2 = pairwise_mean(o2);
        no2 = pairwise_mean(no2);
        co2 = pairwise_mean(co2);
    }

    void init(Logger& logger){
        logger.setAttr("sw", sw);
        logger.setAttr("lw", lw);
    }

    void calculate_radiation(State& state,
                             const std::vector<Superparticle>& superparticles
                             ) {
        if (lw || sw) {
            if (first) {
                prepare_rad_solver_input(state);
                first = false;
                nlay = T_lay_app.size();
            }
            std::vector<double> ch4vmr(nlay, 0);
            std::vector<double> cfc11vmr(nlay, 0);
            std::vector<double> cfc12vmr(nlay, 0);
            std::vector<double> cfc22vmr(nlay, 0);
            std::vector<double> ccl4vmr(nlay, 0);
            std::vector<double> h2o(nlay, 0);
            std::vector<double> o3(nlay, 0);
            std::vector<double> o2(nlay, 0);
            std::vector<double> co2(nlay, 0);
            std::vector<double> no2(nlay, 0);

            std::vector<double> cliqwp(nlay, 0);
            std::vector<double> reliq(nlay, 0);
            std::vector<double> hr(nlay, 0);

            double** uflxlw;
            double** uflxsw;
            double** dflxlw;
            double** dflxsw;
            double** hrlw;
            double** hrsw;

            calculate_cloudproperties(superparticles, state.grid, cliqwp, reliq);

            if (lw) {
                cfpda_rrtm_lw_cld(
                    1, nlay, p_lvl_app.data(), T_lay_app.data(), h2o.data(),
                    o3.data(), co2.data(), ch4vmr.data(), no2.data(), o2.data(),
                    cfc11vmr.data(), cfc12vmr.data(), cfc22vmr.data(),
                    ccl4vmr.data(), cliqwp.data(), reliq.data(), &uflxlw,
                    &dflxlw, &hrlw);
                std::transform(&hrlw[0][0], &hrlw[0][0] + nlay, hr.begin(),
                               hr.begin(), std::plus<void>());
            }
            if (sw) {
                cfpda_rrtm_sw_cld(
                    1, nlay, p_lvl_app.data(), T_lay_app.data(), h2o.data(),
                    o3.data(), co2.data(), ch4vmr.data(), no2.data(), o2.data(),
                    cliqwp.data(), reliq.data(), &uflxsw, &dflxsw, &hrsw);
                std::transform(&hrsw[0][0], &hrsw[0][0] + nlay, hr.begin(),
                               hr.begin(), std::plus<void>());
            }

            std::vector<double> Enet;
            for (unsigned int i = nlay - 1;
                 i > (nlay - 1 - state.layers.size()); --i) {
                Enet.push_back(-hr[i] / (24. * 60. * 60.) * state.grid.length * C_P *
                               RHO_AIR);
            }
            std::copy(Enet.begin(), Enet.end(),
                      member_iterator(state.layers.begin(), &Layer::E));
        }
    }

    void prepare_rad_solver_input(State& state) {
        // append afglus pressure to model pressure
        double p_ref = state.levels.back().p / 100.;
        index = std::lower_bound(p.begin(), p.end(), p_ref) - p.begin();

        std::vector<double> p_buf;
        std::reverse_copy(p.begin(), p.begin() + index,
                          std::back_inserter(p_buf));
        std::transform(p_buf.begin(), p_buf.end(), p_buf.begin(),
                       [](double a) { return a * 100; });

        p_lvl_app =
            concatonate(member_iterator(state.levels.begin(), &Level::p),
                        member_iterator(state.levels.end(), &Level::p),
                        p_buf.begin(), p_buf.end());

        std::transform(p_lvl_app.begin(), p_lvl_app.end(), p_lvl_app.begin(),
                       [](double a) { return a / 100; });
        std::reverse(p_lvl_app.begin(), p_lvl_app.end());
        std::vector<double> T_buf;
        std::reverse_copy(T.begin(), T.end(), std::back_inserter(T_buf));
        T_buf = pairwise_mean(T_buf);
        T_lay_app =
            concatonate(member_iterator(state.layers.begin(), &Layer::T),
                        member_iterator(state.layers.end(), &Layer::T),
                        T_buf.begin() + T_buf.size() - index, T_buf.end());
        std::reverse(T_lay_app.begin(), T_lay_app.end());
    }

    bool sw;
    bool lw;
    std::vector<double> T_lay_app;
    std::vector<double> p_lvl_app;
    std::vector<double> z;
    std::vector<double> p;
    std::vector<double> T;
    std::vector<double> h2o;
    std::vector<double> o3;
    std::vector<double> o2;
    std::vector<double> co2;
    std::vector<double> no2;
    std::vector<double> air;
    int index;
    int nlay;
    bool first = true;
};
