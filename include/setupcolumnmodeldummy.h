#pragma once
#include <random>
#include "columnmodel.h"
#include "setupstate.h"
#include "grid.h"

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource() {
    return mkSPSCH<OIt>(500, 0, 1.e+8,
                        std::lognormal_distribution<double>(log(0.02e-6), 1.5),
                        std::mt19937_64());
}

State createState(Grid grid) {
    double t = 0;
    std::vector<double> z_lay = grid.getlays();

    return {t,
            {std::vector<double>(grid.n_lvl, 5), grid.length},
            {hydrostatic_pressure(z_lay, 100000), grid.length},
            {linear_temperature(z_lay, 288), grid.length},
            {std::vector<double>(grid.n_lay, 0), grid.length},
            {exponential_qv(z_lay, 0.015, 2000), grid.length}};
}

ColumnModel createColumnModel() {
    Grid grid{3000., 100.};
    double t_max = 100;
    double dt = 1;
    double dz = 50;
    return ColumnModel(createState(grid),
                       createParticleSource<ColumnModel::OIt>(), t_max, dt, dz,
                       grid);
}
