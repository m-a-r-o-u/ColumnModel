#include <iostream>
#include <iterator>
#include <random>
#include "columnmodel.h"
#include "thermodynamic.h"
#include "logger.h"
#include "superparticle_population.h"

void ColumnModel::run(std::shared_ptr<Logger> logger) {
    logger->log(state, superparticles.population, grid);
    while (is_running()) {
        step();
        logger->log(state, superparticles.population, grid);
    }
}

void ColumnModel::step() {
    source->generateParticles(std::back_inserter(superparticles.population), dt);

    State old_state(state);
    for (auto& superparticle : superparticles.population) {
        double w = old_state.w(superparticle.z, superparticles.dz);
        double T = old_state.T(superparticle.z, superparticles.dz);
        double p = old_state.p(superparticle.z, superparticles.dz);
        double qv = old_state.qv(superparticle.z, superparticles.dz);
        double E = old_state.E(superparticle.z, superparticles.dz);
        double S = saturation(T, p, qv);

        if (superparticle.is_nucleated ||
            will_nucleate(superparticle.r_dry, S, T)) {
            auto tendencies = calc_tendencies(superparticle, w, S, T, E, dt);
            apply_tendencies_to_superparticle(superparticle, tendencies, w);
            apply_tendencies_to_state(superparticle, tendencies);
        }
    }
    state.qv.advect(state.w, dt);
}

void ColumnModel::apply_tendencies_to_superparticle(
    Superparticle& superparticle, Tendencies& tendencies, const double w) {
    superparticle.qc += tendencies.dqc;
    superparticle.z += w * dt;
    nucleation(superparticle);
}

void ColumnModel::apply_tendencies_to_state(Superparticle& superparticle,
                                            const Tendencies& tendencies) {
    state.T.change(superparticle.z, superparticles.dz, tendencies.dT);
    state.qv.change(superparticle.z, superparticles.dz, -tendencies.dqc);
}

Tendencies ColumnModel::calc_tendencies(const Superparticle& superparticle,
                                        const double w, const double S,
                                        const double T, const double E,
                                        const double dt) {
    auto tendencies = condensation(superparticle.qc, superparticle.N,
                                   superparticle.r_dry, S, T, E, dt);
    return tendencies;
}

void ColumnModel::nucleation(Superparticle& superparticle) {
    auto r = radius(superparticle.qc, superparticle.N, superparticle.r_dry);
    if ((r - superparticle.r_dry) / superparticle.r_dry < 1.e-5) {
        superparticle.is_nucleated = false;
    } else {
        superparticle.is_nucleated = true;
    }
}

bool ColumnModel::is_running() {
    runs++;
    state.t = runs * dt;
    if (state.t < t_max) {
        return true;
    } else {
        return false;
    }
}
