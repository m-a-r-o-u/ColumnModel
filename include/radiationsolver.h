#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "readatm_utils.h"
#include "backgroundlevel.h"
#include "state.h"
#include "fpda_rrtm_lw_cld.h"

template <typename IIt, typename IIt2>
std::vector<double> concatonate(IIt a_first, IIt a_last, IIt2 b_first,
                                IIt2 b_last) {
    std::vector<double> res;
    res.reserve(a_last - a_first + b_last - b_first);
    res.insert(res.end(), a_first, a_last);
    res.insert(res.end(), b_first, b_last);
    return res;
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

    void lw(State& state, Superparticles& superparticles) {
        if (first) {
            prepare_rad_solver_input(state);
            first = false;
            nlay = T_lay_app.size();
        }

        //nlay = T.size();
        std::vector<double> ch4vmr(nlay, 0);
        std::vector<double> cfc11vmr(nlay, 0);
        std::vector<double> cfc12vmr(nlay, 0);
        std::vector<double> cfc22vmr(nlay, 0);
        std::vector<double> ccl4vmr(nlay, 0);

        std::vector<double> cliqwp(nlay, 0);
        std::vector<double> reliq(nlay, 0);
        //std::vector<double> cliqwp_dummy = calculate_cloudliquidwaterpath();
        //std::vector<double> reliq(nlay, 0);

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

        int width = 15;
        std::cout << std::endl;
        std::cout << std::setw(width) << "uflx" << std::setw(width) << "dflx"
                  << std::setw(width) << "hr" << std::endl;
        for (int i = 0; i < nlay + 1; ++i) {
            std::cout << std::setw(width) << uflx[0][i] << std::setw(width)
                      << dflx[0][i];
            if (i >= nlay) {
                std::cout << std::setw(width) << "No Value";
            } else {
                std::cout << std::setw(width) << hr[0][i];
            }
            std::cout << std::endl;
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
        // append afglus temperature to model temperature
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
