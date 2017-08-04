#pragma once
#include <iterator>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
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

    std::vector<double> qc_sum(grid.n_lay, 0);
    std::vector<double> r_2(grid.n_lay, 0);
    std::vector<double> r_3(grid.n_lay, 0);
    std::vector<double> r_eff(grid.n_lay, 0);

    for (auto sp : superparticles) {
        std::vector<double>::iterator ubound =
            std::upper_bound(lvls.begin(), lvls.end(), sp.z);
        assert(ubound != lvls.begin());

        int index = std::distance(lvls.begin(), ubound) - 1;
        assert(index >= 0);
        assert(index < grid.n_lay);
        if (sp.is_nucleated) {
            qc_sum[index] += sp.qc;
            r_2[index] += std::pow(radius(sp.qc, sp.N, sp.r_dry), 2);
            r_3[index] += std::pow(radius(sp.qc, sp.N, sp.r_dry), 3);
        }
    }

    std::transform(r_3.begin(), r_3.end(), r_2.begin(), r_eff.begin(),
                   std::divides<void>());

    std::replace_if(r_eff.begin(), r_eff.end(),
                    [](const double& a) { return std::isnan(a); }, 0);

    std::transform(
        r_eff.begin(), r_eff.end(), r_eff.begin(),
        std::bind(std::multiplies<void>(), std::placeholders::_1, 1.e6));

    std::transform(
        qc_sum.begin(), qc_sum.end(), qc_sum.begin(),
        std::bind(std::multiplies<void>(), std::placeholders::_1, 1.e3));

    std::transform(r_eff.begin(), r_eff.end(), reliq.begin(), reliq.begin(),
                   std::plus<void>());

    std::transform(qc_sum.begin(), qc_sum.end(), cliqwp.begin(), cliqwp.begin(),
                   std::plus<void>());

    std::transform(
        cliqwp.begin(), cliqwp.end(), cliqwp.begin(),
        std::bind(std::multiplies<void>(), std::placeholders::_1, 50.));

    std::reverse(reliq.begin(), reliq.end());
    std::reverse(cliqwp.begin(), cliqwp.end());

//    int width = 12;
//    std::cout << std::endl;
//    std::cout << std::setw(width) << "qc";
//    std::cout << std::setw(width) << "reff";
//    std::cout << std::endl;
//    for (int i = 0; i < grid.n_lay - 1; ++i) {
//        std::cout << std::setw(width) << qc_sum.at(i);
//        std::cout << std::setw(width) << r_eff[i];
//        std::cout << std::endl;
//    }
}

struct RadiationSolver {
   public:
    RadiationSolver(std::string filename) {
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

    void lw(State& state, std::vector<Superparticle>& superparticles,
            Grid grid) {
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

        std::vector<double> cliqwp(nlay, 0);
        std::vector<double> reliq(nlay, 0);
        calculate_cloudproperties(superparticles, grid, cliqwp, reliq);
        double re_min = 2.5;
        double re_max = 60.;
        std::replace_if(reliq.begin(), reliq.end(),
                        [re_min](double a) { return (a > 0 && a < re_min); },
                        re_min);
        std::replace_if(reliq.begin(), reliq.end(),
                        [re_max](double a) { return (a > re_max); },
                        re_max);

        std::vector<double> h2o_dummy(nlay, 0);
        std::vector<double> o3_dummy(nlay, 0);
        std::vector<double> o2_dummy(nlay, 0);
        std::vector<double> co2_dummy(nlay, 0);
        std::vector<double> no2_dummy(nlay, 0);

        double** uflx;
        double** dflx;
        double** hr;

        cfpda_rrtm_lw_cld(
            1, nlay, p_lvl_app.data(), T_lay_app.data(), h2o_dummy.data(),
            o3_dummy.data(), co2_dummy.data(), ch4vmr.data(), no2_dummy.data(),
            o2_dummy.data(), cfc11vmr.data(), cfc12vmr.data(), cfc22vmr.data(),
            ccl4vmr.data(), cliqwp.data(), reliq.data(), &uflx, &dflx, &hr);

        std::vector<double> Enet;

        for (unsigned int i = nlay - 1; i > (nlay - 1 - state.layers.size());
             --i) {
            Enet.push_back(- hr[0][i] / (24. * 60. * 60.) * grid.length * C_P *
                           RHO_AIR);
        }
        std::copy(Enet.begin(), Enet.end(),
                  member_iterator(state.layers.begin(), &Layer::E));

        //        int width = 15;
        //        std::cout << std::endl;
        //        std::cout << std::setw(width) << "r_eff";
        //        std::cout << std::setw(width) << "cliqwp";
        //        std::cout << std::setw(width) << "p_lvl_app";
        //        std::cout << std::endl;
        //        for (int i = 0; i < nlay + 1; ++i) {
        //            std::cout << std::setw(width) << reliq[i];
        //            std::cout << std::setw(width) << cliqwp[i];
        //            std::cout << std::setw(width) << p_lvl_app[i];
        //            std::cout << std::endl;
        //        }
        //
        //        std::cout << std::endl;
        //        std::cout << std::setw(width) << "uflx" << std::setw(width) <<
        //        "dflx"
        //                  << std::setw(width) << "hr" << std::endl;
        //        for (int i = 0; i < nlay + 1; ++i) {
        //            std::cout << std::setw(width) << uflx[0][i] <<
        //            std::setw(width)
        //                      << dflx[0][i];
        //            if (i >= nlay) {
        //                std::cout << std::setw(width) << "No Value";
        //            } else {
        //                std::cout << std::setw(width) << hr[0][i];
        //            }
        //            std::cout << std::endl;
        //        }
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
