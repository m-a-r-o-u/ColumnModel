#include "columnmodel.h"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <random>
#include "analize_sp.h"
#include "collision.h"
#include "continuouse_state_view.h"
#include "logger.h"
#include "member_iterator.h"
#include "saturation_fluctuations.h"
#include "thermodynamic.h"
#include "twomey.h"

void cooling_the_column(State& state, double dt) {
    std::vector<double> cooling(state.grid.n_lay, -2.3e-5 * 24. * dt);
    std::transform(member_iterator(state.layers.begin(), &Layer::T),
                   member_iterator(state.layers.end(), &Layer::T),
                   cooling.begin(),
                   member_iterator(state.layers.begin(), &Layer::T),
                   [](double x, double y) { return x + y; });
}

void check_state(State& state) {
    for (auto l : state.layers) {
        if (l.qv < 0) {
            std::cout << "index qv" << std::endl;
            for (unsigned int i = 0; i < state.layers.size(); ++i) {
                std::cout << i << " " << state.layers[i].qv << std::endl;
            }
            throw std::logic_error("qv is smaller then zero ");
        }
    }
}

void check_S(const State& s, const State& os) {
    for (unsigned int i = 0; i > s.layers.size(); ++i) {
        Layer lay = s.layers[i];
        Layer olay = os.layers[i];
        double S = super_saturation(lay.T, lay.p, lay.qv);
        double oS = super_saturation(olay.T, olay.p, olay.qv);
        if (oS > 0) {
            if (S < 0) {
                throw std::logic_error(
                    "super saturation is smalle zero, this is not correct if "
                    "radiation solver and fluctuation solver a turned off");
            }
        }
    }
}

void check_sp(const Superparticle& sp) {
    if (sp.qc < 0.) {
        throw std::logic_error("qc of one superparticle is smaller zero: " +
                               std::to_string(sp.qc));
    }
    if (sp.N < 0.) {
        throw std::logic_error("N of one superparticle is smaller zero: " +
                               std::to_string(sp.N));
    }
}

void check_superparticles(const std::vector<Superparticle>& sp,
                          const Grid& grid) {
    for (auto s : sp) {
        check_sp(s);
    }
}

void ColumnModel::run(std::shared_ptr<Logger> logger) {
    logger->initialize(state, dt);
    radiation_solver.init(*logger);
    source->init(*logger);

    logger->log(state, superparticles);
    while (is_running()) {
        step();
        log_every_seconds(logger, 30.);
    }
}

void ColumnModel::step() {
    advection_solver->advect(state, dt);
    advection_solver->setupdraft(state, runs * dt);
    advection_solver->keepcloudbase(state);

    fluctuations->refresh(superparticles);

    State old_state(state);

    source->generateParticles(std::back_inserter(superparticles), state, dt,
                              superparticles);

    if (true) {
        check_state(state);
        check_superparticles(superparticles, state.grid);
    }
    do_condensation(old_state);
    removeUnnucleated(superparticles);
    for (const auto& sp : superparticles) {
        assert(sp.is_nucleated == true);
        assert(sp.qc >= 0);
        assert(sp.z >= 0);
        assert(sp.N >= 0);
    }
    do_collisions();
    removeUnnucleated(superparticles);
    radiation_solver.calculate_radiation(state, superparticles);
}

void ColumnModel::do_condensation(State& old_state) {
    for (auto& sp : superparticles) {
        Layer lay = old_state.layer_at(sp.z);
        Level lvl = old_state.upper_level_at(sp.z);
        double S = super_saturation(lay.T, lay.p, lay.qv) +
                   fluctuations->getFluctuation(sp, dt);
        if (sp.is_nucleated) {
            auto tendencies = calc_tendencies(sp, S, lay.T, lay.E, dt);
            apply_tendencies_to_superparticle(sp, tendencies, lvl);
            apply_tendencies_to_state(sp, tendencies);
        }
    }
}

void ColumnModel::do_collisions() {
    if (collisions->needs_sorted_superparticles()) {
        std::sort(superparticles.begin(), superparticles.end(),
                  [](const auto& a, const auto& b) { return a.z < b.z; });
    }
    auto collision_tendencies =
        collisions->collide(superparticles, state.grid, dt);
    for (const auto& c : collision_tendencies) {
        if (std::isnan(c.dqc)) {
            throw std::logic_error("collison dqc is nan");
        }
    }
    apply_collision_tendencies(superparticles, collision_tendencies);
}

void ColumnModel::apply_collision_tendencies(
    std::vector<Superparticle>& sps,
    const std::vector<SpMassTendencies>& tendencies) {
    assert(sps.size() == tendencies.size());
    for (size_t i = 0; i < sps.size(); ++i) {
        sps[i].N += tendencies[i].dN;
        sps[i].qc += tendencies[i].dqc;
        sps[i].update();
    }
}

void ColumnModel::log_every_seconds(std::shared_ptr<Logger> logger,
                                    double dt_out) {
    if (!std::abs(std::remainder(runs * dt, dt_out))) {
        logger->log(state, superparticles);
    }
}

void ColumnModel::apply_tendencies_to_superparticle(Superparticle& sp,
                                                    Tendencies& tendencies,
                                                    const Level& lvl) {
    sp.v = lvl.w - sedimentation->fall_speed(sp.radius());
    double cfl = sp.v * dt / state.grid.length;
    if (cfl > 1) {
        throw std::logic_error("the cfl criteria is broken: cfl=" +
                               std::to_string(cfl));
    }
    sp.z += dt * sp.v;
    sp.qc += tendencies.dqc;
    sp.update();
    check_sp(sp);
}

void ColumnModel::apply_tendencies_to_state(const Superparticle& sp,
                                            const Tendencies& tendencies) {
    if (sp.z >= 0) {
        state.change_layer(sp.z, {0, 0, -tendencies.dqc, 0});
    }
    if (sp.z <= 0. && sp.qc > 0 && sp.N > 0) {
        state.qr_ground += sp.qc;
    }
}

Tendencies ColumnModel::calc_tendencies(const Superparticle& superparticle,
                                        const double S, const double T,
                                        const double E, const double dt) {
    auto tendencies = condensation(superparticle.qc, superparticle.N,
                                   superparticle.r_dry, S, T, E, dt);
    return tendencies;
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
