#pragma once
#include <yaml-cpp/yaml.h>
#include <memory>
#include <random>
#include "advect.h"
#include "columnmodel.h"
#include "grid.h"
#include "layer_quantities.h"
#include "level_quantities.h"
#include "radiationsolver.h"
#include "setupstate.h"
#include "state.h"
#include "twomey.h"

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(
    const Grid& grid, const YAML::Node& config) {

    int N_sp = config["N_sp"].as<int>();
    int N_lay = grid.n_lay;
    return mkTwomey<OIt>(N_sp, N_lay);
}

//std::function<State(const Grid&)> createState(const YAML::Node& config) {
//    return [](const Grid& grid) { return State(...); }
//}

State createState(const Grid& grid, const YAML::Node& config) {
    double T0 = config["T0"].as<double>();
    double p0 = config["p0"].as<double>();
    double cloud_base = config["cloud_base"].as<double>();
    double w = config["w"].as<double>();

    State state{0, {}, {}, grid, cloud_base, w};

    for (const auto& el : grid.getlays()) {
        state.layers.push_back(
            {linear_temperature(el, T0), hydrostatic_pressure(el, p0), 0, 0});
    }

    int index = std::floor(cloud_base / grid.length) - 1;
    state.layers[index].qv =
        saturation_vapor(state.layers[index].T, state.layers[index].p);

    for (const auto& el : grid.getlvls()) {
        state.levels.push_back({w, hydrostatic_pressure(el, p0)});
    }
    return state;
}

RadiationSolver createRadiationSolver(const YAML::Node& config) {
    bool sw = config["sw"].as<bool>();
    bool lw = config["lw"].as<bool>();
    std::string data_path = config["data_path"].as<std::string>();
    return RadiationSolver(data_path, sw, lw);
}

std::unique_ptr<Advect> createAdvectionSolver(const YAML::Node& config) {
    return mkFirstOrder();
}

std::unique_ptr<Grid> createGrid(const YAML::Node& config) {
    double toa = config["toa"].as<double>();
    double gridlength = config["gridlength"].as<double>();
    return std::make_unique<Grid>(toa, gridlength);
}

ColumnModel createColumnModel(const YAML::Node& config) {
    double t_max = config["t_max"].as<double>();
    double dt = config["dt"].as<double>();

    auto grid = createGrid(config["grid"]);

    auto advection_solver = createAdvectionSolver(config["advection"]);
    auto state = createState(*grid, config["initial_state"]);
    auto radiation_solver = createRadiationSolver(config["radiation"]);
    auto source = createParticleSource<ColumnModel::OIt>(
        *grid, config["particle_source"]);

    return ColumnModel(state, std::move(source), t_max, dt, radiation_solver,
                       std::move(grid), std::move(advection_solver));
}
