#pragma once
#include <random>
#include "columnmodel.h"
#include "grid.h"
#include "layer_quantities.h"
#include "level_quantities.h"
#include "radiationsolver.h"
#include "setupstate.h"
#include "state.h"
#include "twomey.h"

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(double z_insert,
                                                               int N,
                                                               int N_multi,int N_lay) {
    if (false) {
        return mkSPSCH<OIt>(
            z_insert, N, N_multi,
            std::lognormal_distribution<double>(log(0.02e-6), 1.5),
            std::mt19937_64());
    }
    if (true) {
        return mkTwomey<OIt>(N, N_multi, N_lay);
    }
}

State createState(Grid grid, double w, double p0) {
    State state{0, {}, {}, grid};

    for (const auto& el : grid.getlays()) {
        state.layers.push_back(
            {linear_temperature(el, 286), hydrostatic_pressure(el, p0), 0, 0});
    }

    int index = std::floor(grid.z0 / grid.length) - 1;
    state.layers[index].qv =
        saturation_vapor(state.layers[index].T, state.layers[index].p);
    //state.layers[index].qv += 0.001;

    for (const auto& el : grid.getlvls()) {
        state.levels.push_back({w, hydrostatic_pressure(el, p0)});
    }
    return state;
}

RadiationSolver createRadiationSolver() {
    return RadiationSolver(
        "/home/m/Mares.Barekzai/phd/projects/column_model/data/afglus.dat");
}

ColumnModel createColumnModel() {
    double t_max = 30 * 60;
    double dt = 0.1;
    double w = 1;
    double gridlength = 100.;
    double z0 = 1000;
    double toa = z0 + w * t_max + gridlength*5;
    int N = 1.e6;
    assert(N >= 1);
    int N_multi = 1.e8 / double(N);
    double p0 = 100000;
    Grid grid{toa, gridlength, z0};

    bool lw = false;
    bool sw = false;

    std::cout << grid.n_lay << std::endl;
    auto state = createState(grid, w, p0);
    auto radiation_solver = createRadiationSolver();

    return ColumnModel(
        state, createParticleSource<ColumnModel::OIt> (z0, N, N_multi, grid.n_lay),
        t_max, dt, radiation_solver, lw, sw, N);
}
