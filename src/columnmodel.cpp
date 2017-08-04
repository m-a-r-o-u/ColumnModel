#include "columnmodel.h"
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <random>
#include "continuouse_state_view.h"
#include "logger.h"
#include "member_iterator.h"
#include "thermodynamic.h"

void ColumnModel::log_every_seconds(std::shared_ptr<Logger> logger,
                                    double dt_out) {
    if (!std::abs(std::remainder(runs * dt, dt_out))) {
        logger->log(state, superparticles, grid);
    }
}

void ColumnModel::run(std::shared_ptr<Logger> logger) {
    logger->log(state, superparticles, grid);
    while (is_running()) {
        step();
        log_every_seconds(logger, 60.);
    }
}

void ColumnModel::step() {
    source->generateParticles(std::back_inserter(superparticles), dt);
    State old_state(state);

    for (auto& sp : superparticles) {
        Layer lay = old_state.layer_at(sp.z, grid);
        Level lvl = old_state.upper_level_at(sp.z, grid);

        double S = saturation(lay.T, lay.p, lay.qv);

        if (sp.is_nucleated || will_nucleate(sp.r_dry, S, lay.T)) {
            auto tendencies = calc_tendencies(sp, S, lay.T, lay.E, dt);
            apply_tendencies_to_superparticle(sp, tendencies, lvl.w);
            apply_tendencies_to_state(sp, tendencies);
        }
        sp.z += lvl.w * dt - dt * fall_speed(radius(sp.qc, sp.N, sp.r_dry));
    }

    auto cloud_bottom_lay =
        state.layers.begin() + std::floor(grid.z0 / grid.length) - 1;
    auto cloud_bottom_lev =
        state.levels.begin() + std::floor(grid.z0 / grid.length) - 1;
    advect_first_order(member_iterator(cloud_bottom_lay, &Layer::qv),
                       member_iterator(state.layers.end(), &Layer::qv),
                       member_iterator(cloud_bottom_lev, &Level::w),
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
    if (state.t < t_max) {
        return true;
    } else {
        return false;
    }
}
