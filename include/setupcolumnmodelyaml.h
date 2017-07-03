#pragma once
#include "yaml-cpp/yaml.h"
#include <random>
#include <string>
#include <vector>
#include "columnmodel.h"
#include <stdexcept>
#include "linearfield.h"
#include "constants.h"

void throw_if_config_argument_invalid(const YAML::Node& config, std::string s) {
    if (!config[s]) {
        throw std::invalid_argument(s + " does not specify a YAML::Node");
    }
}

template <typename OIt>
std::unique_ptr<SuperParticleSource<OIt>> createParticleSource(
    const YAML::Node& config) {
    return mkSPSCH<OIt>(500, 1, 1.e+8,
                        std::lognormal_distribution<double>(log(0.02e-6), 1.5),
                        std::mt19937_64());
}

std::vector<double> exponential_qv(std::vector<double> z, double qv0,
                                   double zc) {
    std::vector<double> qv;
    for (auto el : z) {
        qv.push_back(qv0 * std::exp(-el / zc));
    }
    return qv;
}

std::vector<double> linear_temperature(const std::vector<double>& z,
                                       double T0) {
    std::vector<double> T;
    for (auto el : z) {
        T.push_back(T0 - LAPSE_RATE_A * el);
    }
    return T;
}
std::vector<double> hydrostatic_pressure(const std::vector<double>& z,
                                         double p0) {
    std::vector<double> p;
    for (auto el : z) {
        p.push_back(p0 - G * el);
    }
    return p;
}

std::vector<double> create_z(double dz, double n) {
    std::vector<double> z;
    for (int i = 0; i < n; ++i) {
        z.push_back(i * dz);
    }
    return z;
}

State createState(const YAML::Node& config) {
    throw_if_config_argument_invalid(config, "initial_conditions");
    const YAML::Node initial_conditions = config["initial_conditions"];
    throw_if_config_argument_invalid(initial_conditions, "initial_state");
    const YAML::Node initial_state = initial_conditions["initial_state"];
    throw_if_config_argument_invalid(initial_state, "w");
    const YAML::Node w = initial_state["w"];
    const double wind_max = w["max"].as<double>();
    const double total_height = initial_conditions["total_height"].as<double>();
    const double grid_length = initial_conditions["grid_length"].as<double>();
    size_t grid_points = std::ceil(total_height / grid_length);

    const YAML::Node p = initial_state["p"];
    const double p0 = p["p0"].as<double>();

    const YAML::Node E = initial_state["E"];
    const double E0 = E["E0"].as<double>();

    const YAML::Node qv = initial_state["qv"];
    const double qv0 = qv["qv0"].as<double>();
    const double zc = qv["zc"].as<double>();

    std::vector<double> z = create_z(grid_length, grid_points);

    std::vector<double> pressure_values(hydrostatic_pressure(z, p0));
    std::vector<double> wind_values(grid_points, wind_max);
    std::vector<double> radiation_values(grid_points, E0);
    std::vector<double> T_values(linear_temperature(z, T0));
    std::vector<double> qv_values(exponential_qv(z, qv0, zc));

    return {{wind_values, grid_length},
            {pressure_values, grid_length},
            {T_values, grid_length},
            {radiation_values, grid_length},
            {qv_values, grid_length}};
}

ColumnModel createColumnModel(const YAML::Node& config) {
    return ColumnModel(createState(config["initial_conditions"]),
                       createParticleSource<ColumnModel::OIt>(config), 2, 1,
                       50);
}
