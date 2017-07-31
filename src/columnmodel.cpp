#include <iostream>
#include <iterator>
#include <random>
#include "columnmodel.h"
#include "thermodynamic.h"
#include "logger.h"
#include "continuouse_state_view.h"
#include "member_iterator.h"
#include <cstdlib>

void ColumnModel::run(std::shared_ptr<Logger> logger) {
    logger->log(state, superparticles, grid);
    while (is_running()) {
        step();
        logger->log(state, superparticles, grid);
    }
}

void ColumnModel::step() {
    source->generateParticles(std::back_inserter(superparticles), dt);
    State old_state(state);

    for (auto& superparticle : superparticles) {
        Layer lay = old_state.layer_at(superparticle.z, grid);
        Level lvl = old_state.upper_level_at(superparticle.z, grid);

        double S = saturation(lay.T, lay.p, lay.qv);

        if (superparticle.is_nucleated ||
            will_nucleate(superparticle.r_dry, S, lay.T)) {
            auto tendencies =
                calc_tendencies(superparticle, S, lay.T, lay.E, dt);
            apply_tendencies_to_superparticle(superparticle, tendencies, lvl.w);
            apply_tendencies_to_state(superparticle, tendencies);
        }
        superparticle.z +=
            lvl.w * dt +
            dt * fall_speed(radius(superparticle.qc, superparticle.N,
                                   superparticle.r_dry));
    }
    advect_first_order(member_iterator(state.layers.begin(), &Layer::qv),
                       member_iterator(state.layers.end(), &Layer::qv),
                       member_iterator(state.levels.begin(), &Level::w),
                       grid.length, dt);
    radiation_solver.lw(state, superparticles, grid);
}

void ColumnModel::apply_tendencies_to_superparticle(
    Superparticle& superparticle, Tendencies& tendencies, const double w) {
    superparticle.qc += tendencies.dqc;
    nucleation(superparticle);
}

void ColumnModel::apply_tendencies_to_state(Superparticle& superparticle,
                                            const Tendencies& tendencies) {
    state.change_layer(superparticle.z, grid,
                       //{tendencies.dT, 0, -tendencies.dqc, 0});
                       {0, 0, -tendencies.dqc, 0});
}

Tendencies ColumnModel::calc_tendencies(const Superparticle& superparticle,
                                        const double S, const double T,
                                        const double E, const double dt) {
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
    // if (state.t < t_max) {
    if (runs < 50 * 60) {
        return true;
    } else {
        return false;
    }
}
