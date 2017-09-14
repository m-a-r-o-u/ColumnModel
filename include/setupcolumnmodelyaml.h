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
#include "saturation_fluctuations.h"
#include "setupstate.h"
#include "state.h"
#include "twomey.h"

template <typename OIt, typename G>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(
    G& gen, const Grid& grid, const YAML::Node& config) {
    std::string type = config["type"].as<std::string>();
    int N_sp = config["N_sp"].as<int>();
    int N_lay = grid.n_lay;
    return mkPS<OIt>(gen, N_sp, N_lay, type);
}

void setQvProfile(double base, double roof, State& state) {
    int b_index = std::floor(base / state.grid.length) - 1;
    int r_index = std::floor(roof / state.grid.length) - 1;
    for (int i = b_index; i <= r_index; ++i){
        state.layers[i].qv =
           saturation_vapor(state.layers[i].T, state.layers[i].p);
    }
}

State createState(const Grid& grid, const YAML::Node& config) {
    double T0 = config["T0"].as<double>();
    double p0 = config["p0"].as<double>();
    double cloud_base = config["cloud_base"].as<double>();
    double cloud_roof = config["cloud_roof"].as<double>();
    double w = config["w"].as<double>();

    State state{0, {}, {}, grid, cloud_base, w};

    for (const auto& el : grid.getlays()) {
        state.layers.push_back(
            {linear_temperature(el, T0), hydrostatic_pressure(el, p0), 0, 0});
    }

    setQvProfile(cloud_base, cloud_roof, state);
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

template <typename G>
std::unique_ptr<FluctuationSolver> createFluctuationSolver(
    G& gen, const Grid& grid, const YAML::Node& config) {
    std::string type = config["type"].as<std::string>();
    double epsilon = config["epsilon"].as<double>();
    return mkFS(gen, type, epsilon, grid);
}

template <typename G>
ColumnModel createColumnModel(G& gen, const YAML::Node& config) {
    double t_max = config["t_max"].as<double>();
    double dt = config["dt"].as<double>();

    auto grid = createGrid(config["grid"]);

    auto advection_solver = createAdvectionSolver(config["advection"]);
    auto state = createState(*grid, config["initial_state"]);
    auto radiation_solver = createRadiationSolver(config["radiation"]);
    auto source = createParticleSource<ColumnModel::OIt>(
        gen, *grid, config["particle_source"]);
    auto fluctuations =
        createFluctuationSolver(gen, *grid, config["fluctuations"]);

    return ColumnModel(state, std::move(source), t_max, dt, radiation_solver,
                       std::move(grid), std::move(advection_solver),
                       std::move(fluctuations));
}
