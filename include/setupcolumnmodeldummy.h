#pragma once
#include <random>
#include "columnmodel.h"
#include "setupstate.h"
#include "grid.h"
#include "layer_quantities.h"
#include "level_quantities.h"

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(double z_insert, int N) {
    return mkSPSCH<OIt>(z_insert, N, 1.e+8,
                        std::lognormal_distribution<double>(log(0.02e-6), 1.5),
                        std::mt19937_64());
}

State createState(Grid grid, double w) {
    State state;
    state.t = 0;

    for (const auto& el : grid.getlays()) {
        state.layers.push_back({linear_temperature(el, 288),
                                hydrostatic_pressure(el, 100000),
                                0, 0});
    }
    state.layers[0].qv = 0.010;
    for (const auto& el : grid.getlvls()) {
        state.levels.push_back({w});
    }
    return state;
}

ColumnModel createColumnModel() {
    //double t_max = 60*20;
    double t_max = 60;
    double dt = 0.1;
    double dz = 50;
    double w = 1;
    double gridlength = 50.;
    double toa = w * t_max + (5 + gridlength);
    Grid grid{toa, gridlength};
    int N = w * dt * 20.;

    return ColumnModel(createState(grid, w),
                       createParticleSource<ColumnModel::OIt>(grid.length, N), t_max, dt,
                       dz, grid);
}
