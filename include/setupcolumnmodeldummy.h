#pragma once
#include <random>
#include "columnmodel.h"
#include "setupstate.h"
#include "grid.h"
#include "state.h"
#include "layer_quantities.h"
#include "level_quantities.h"
#include "radiationsolver.h"
#include "fpda_rrtm_lw_cld.h"

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(double z_insert,
                                                               int N) {
    return mkSPSCH<OIt>(z_insert, N, 1.e+8,
                        std::lognormal_distribution<double>(log(0.02e-6), 1.5),
                        std::mt19937_64());
}

State createState(Grid grid, double w, double p0) {
    State state;
    state.t = 0;

    for (const auto& el : grid.getlays()) {
        state.layers.push_back(
            {linear_temperature(el, 288), hydrostatic_pressure(el, p0), 0, 0});
    }
    state.layers[0].qv = 0.010;

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
    //double t_max = 60*20;
    double t_max = 60;
    double dt = 0.1;
    double dz = 50;
    double w = 1;
    double gridlength = 50.;
    double toa = w * t_max + (gridlength + gridlength / 5);
    Grid grid{toa, gridlength};
    int N = w * dt * 20.;
    double p0 = 100000;

    auto radsolver = createRadiationSolver();
    auto state = createState(grid, w, p0);
    radsolver.lw(state);

    return ColumnModel(state,
                       createParticleSource<ColumnModel::OIt>(grid.length, N),
                       t_max, dt, dz, grid);
}
