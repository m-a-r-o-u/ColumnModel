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
#include "collision.h"
#include "sedimentation.h"
#include "constants.h"
#include <exception>

template <typename OIt, typename G>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(
    G& gen, const Grid& grid, const YAML::Node& config) {
    std::string type = config["type"].as<std::string>();
    int N_sp = config["N_sp"].as<int>();
    unsigned int N_lay = grid.n_lay;
    if (type == "twomey") {
        return std::make_unique<Twomey<OIt, G>>(gen, N_sp, N_lay);
    } 
    else if (type == "no") {
        return std::make_unique<NoParticleSource<OIt, G>>(gen, 0., N_lay);
    }
    else{
        throw std::logic_error("the type of the particle source: " + type + " is not found");
    }
}

void setQvProfile(double base, double roof, State& state) {
    int b_index = std::floor(base / state.grid.length) - 1;
    int r_index = std::floor(roof / state.grid.length) - 1;

    for (int i = 0; i < b_index; ++i){
        state.layers[i].qv =
           saturation_vapor(state.layers[i].T, state.layers[i].p) * 0.95;
    }
    for (unsigned int i = b_index; i < state.grid.n_lay; ++i){
        state.layers[i].qv =
           saturation_vapor(state.layers[i].T, state.layers[i].p) * 0.2;
    }

    for (int i = b_index; i <= r_index; ++i){
        state.layers[i].qv =
           saturation_vapor(state.layers[i].T, state.layers[i].p);
    }
}

State createState(const Grid& grid, const YAML::Node& config) {
    double Theta0 = config["Theta0"].as<double>();
    double ALR = config["ALR"].as<double>();
    double p0 = config["p0"].as<double>();
    double cloud_base = config["cloud_base"].as<double>();
    double cloud_roof = config["cloud_roof"].as<double>();
    double w = config["w"].as<double>();

    State state{0, {}, {}, grid, cloud_base, w};

    double Tb = Theta0 * std::pow(hydrostatic_pressure(cloud_base, p0)/p0, R_G/C_P);

    for (const auto& el : grid.getlays()) {
        if(el<cloud_base){
        double p = hydrostatic_pressure(el, p0); 
        double T = Theta0 * std::pow((p/p0), R_G/C_P);
        state.layers.push_back(
            {T, p, 0, 0});
        }
        else{
        state.layers.push_back(
            {Tb - ALR * (el - cloud_base), hydrostatic_pressure(el, p0), 0, 0});
        }
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
    std::string type = config["type"].as<std::string>();
    if (type == "advectandset"){
        double lifetime = config["lifetime"].as<double>();
        return mkASFO(lifetime);
    }
    else if(type == "firstorder"){
        return mkAFO();
    }
    else if(type == "firstorderupwind"){
        double lifetime = config["lifetime"].as<double>();
        return mkAFOU(lifetime);
    }
    else if(type == "secondorderupwind"){
        double lifetime = config["lifetime"].as<double>();
        return mkASOU(lifetime);
    }
    else if(type == "secondfirstorderupwind"){
        double lifetime = config["lifetime"].as<double>();
        return mkASFOU(lifetime);
    }
    else if(type == "thirdorderupwind"){
        double lifetime = config["lifetime"].as<double>();
        return mkATOU(lifetime);
    }
    else if(type == "sixthorderwickerskamarock"){
        double lifetime = config["lifetime"].as<double>();
        return mkASOWK(lifetime);
    }
    else{
        throw std::logic_error("the type of the advection solver: " + type + " is not found");
    }
}

std::unique_ptr<Grid> createGrid(const YAML::Node& config) {
    double toa = config["toa"].as<double>();
    double gridlength = config["gridlength"].as<double>();
    return std::make_unique<Grid>(toa, gridlength);
}

template <typename G>
std::unique_ptr<FluctuationSolver> createFluctuationSolver(
    G& gen, const YAML::Node& config, const Grid& grid) {
    std::string type = config["type"].as<std::string>();
    double epsilon = config["epsilon"].as<double>();
    double l = config["l"].as<double>();
    return mkFS(gen, type, epsilon, l, grid);
}

std::unique_ptr<Collisions> createCollisionSolver(const Sedimentation& sedi,const YAML::Node& config){
    std::string type = config["type"].as<std::string>();
    if ( type == "hall"){
        return mkHCS(sedi);
    }
    else if (type == "no")
    {
        return  mkNCS();
    }
    else{
        throw std::logic_error("the type of the collision solver: " + type + " is not found");
    }
}

std::unique_ptr<Sedimentation> createSedimentationSolver(const YAML::Node& config){
    std::string type = config["type"].as<std::string>();
    if ( type == "lookup"){
        return mkFSLU();
    }
    else if (type == "no")
    {
        return  mkNFS();
    }
    else{
        throw std::logic_error("the type of the sedimentation solver: " + type + " is not found");
    }
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
        createFluctuationSolver(gen, config["fluctuations"], *grid);
    auto sedimentation = createSedimentationSolver(config["sedimentation"]);
    auto collision_solver = createCollisionSolver(*sedimentation, config["collisions"]);

    return ColumnModel(state, std::move(source), t_max, dt, radiation_solver,
                       std::move(grid), std::move(advection_solver),
                       std::move(fluctuations), std::move(collision_solver),
                       std::move(sedimentation));
}
