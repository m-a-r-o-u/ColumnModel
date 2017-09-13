#pragma once
#include <memory>
#include <random>
#include "columnmodel.h"
#include "grid.h"
#include "layer_quantities.h"
#include "level_quantities.h"
#include "radiationsolver.h"
#include "setupstate.h"
#include "state.h"
#include "twomey.h"
#include "advect.h"

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(
                                                               int N_sp,
                                                               int N_sp_lay) {
    if (true) {
        return mkTwomey<OIt>(N_sp, N_sp_lay);
    }
}

State createState(Grid& grid, double w, double p0, int cloud_base) {
    State state{0, {}, {}, grid};

    for (const auto& el : grid.getlays()) {
        state.layers.push_back(
            {linear_temperature(el, 286), hydrostatic_pressure(el, p0), 0, 0});
    }

    int index = std::floor(cloud_base / grid.length) - 1;
    state.layers[index].qv =
        saturation_vapor(state.layers[index].T, state.layers[index].p);
    // state.layers[index].qv += 0.001;

    for (const auto& el : grid.getlvls()) {
        state.levels.push_back({w, hydrostatic_pressure(el, p0)});
    }
    return state;
}

RadiationSolver createRadiationSolver(bool sw, bool lw) {
    return RadiationSolver(
        "/home/m/Mares.Barekzai/phd/projects/column_model/data/afglus.dat", sw,
        lw);
}

std::unique_ptr<Advect> createAdvectionSolver(const double& cloud_base, const double& w_init){
    if(true) {
        return mkFirstOrder(cloud_base, w_init);
    }
}

ColumnModel createColumnModel() {
    double t_max = 30* 60;
    double dt = 0.1;
    double w = 1;
    double gridlength = 10.;
    double cloud_base = 1000;
    double toa = cloud_base + w * t_max + gridlength * 5;
    int N_sp = 1.e3;
    double p0 = 100000;
    auto grid = std::make_unique<Grid>(toa, gridlength);

    bool lw = false;
    bool sw = false;

    auto advection_solver = createAdvectionSolver(cloud_base, w);
    auto state = createState(*grid, w, p0, cloud_base);
    auto radiation_solver = createRadiationSolver(sw, lw);
    auto source =
        createParticleSource<ColumnModel::OIt>(N_sp, grid->n_lay);

    return ColumnModel(state, std::move(source), t_max, dt, radiation_solver,
                       std::move(grid), std::move(advection_solver));
}
