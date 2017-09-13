#include "columnmodel.h"
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <random>
#include "continuouse_state_view.h"
#include "logger.h"
#include "member_iterator.h"
#include "thermodynamic.h"
#include "twomey.h"
#include "analize_sp.h"
#include "saturation_fluctuations.h"

void cooling_the_column(State& state, double dt){
    std::vector<double> cooling(state.grid.n_lay, -.5e-2 * 1 * dt);
    std::transform(member_iterator(state.layers.begin(), &Layer::T),
              member_iterator(state.layers.end(), &Layer::T),
              cooling.begin(),
              member_iterator(state.layers.begin(), &Layer::T),
              [](double x, double y){return x + y;});
}

void check_state(State& state) {
    for (auto l : state.layers) {
        if (l.qv < 0) {
            std::cout << "qv is smaller then zero: " << std::endl;
            std::copy(member_iterator(state.layers.begin(), &Layer::qv),
                      member_iterator(state.layers.end(), &Layer::qv),
                      std::ostream_iterator<double>(std::cout, "\n"));
            std::exit(0);
        }
    }
}

void check_S(const State& s, const State& os){
    for (unsigned int i = 0; i > s.layers.size(); ++i) {
        Layer lay = s.layers[i];
        Layer olay = os.layers[i];
        double S = super_saturation(lay.T, lay.p, lay.qv);
        double oS = super_saturation(olay.T, olay.p, olay.qv);
        if ( oS > 0) {
            if ( S < 0){
                std::cout << "The saturation condenses smaller then zero" << std::endl;
                std::exit(0);
            }
        }
    }
}

void check_superparticles(std::vector<Superparticle>& sp) {
    for (auto s : sp) {
        if (s.N > sp[0].N || s.N < sp[0].N) {
            std::exit(0);
        }
        if (s.qc < 0.) {
            std::cout << s << std::endl;
            std::cout << "qc is smaller then 0 problemo" << std::endl;
            std::exit(0);
        }
    }
}

void ColumnModel::run(std::shared_ptr<Logger> logger) {
    logger->initialize(state, dt);
    radiation_solver.init(*logger);
    source->init(*logger);
    advection_solver->init(*logger);

    logger->log(state, superparticles);
    while (is_running()) {
        step();
        log_every_seconds(logger, 60.);
    }
}

void ColumnModel::step() {
    advection_solver->advect(state, dt);
    source->generateParticles(std::back_inserter(superparticles), state, dt,
                              superparticles);
    fluctuations->refresh(superparticles);

    State old_state(state);

    if (true) {
        check_state(state);
        check_superparticles(superparticles);
    }

    for (auto& sp : superparticles) {
        Layer lay = old_state.layer_at(sp.z);
        Level lvl = old_state.upper_level_at(sp.z);

        double S = super_saturation(lay.T, lay.p, lay.qv) + fluctuations->getFluctuation(sp, dt);

        if (sp.is_nucleated) {
            auto tendencies = calc_tendencies(sp, S, lay.T, lay.E, dt);
            apply_tendencies_to_superparticle(sp, tendencies);
            apply_tendencies_to_state(sp, tendencies);
            sp.z += lvl.w * dt - dt * fall_speed(radius(sp.qc, sp.N, sp.r_dry));
        }
    }

    radiation_solver.calculate_radiation(state, superparticles);

    removeUnnucleated(superparticles);
    if (true){
        check_S(state, old_state);
    }
}

void ColumnModel::log_every_seconds(std::shared_ptr<Logger> logger,
                                    double dt_out) {
    if (!std::abs(std::remainder(runs * dt, dt_out))) {
        logger->log(state, superparticles);
    }
}


void ColumnModel::apply_tendencies_to_superparticle(
    Superparticle& sp, Tendencies& tendencies) {
    sp.qc += tendencies.dqc;
    nucleation(sp);
}

void ColumnModel::apply_tendencies_to_state(Superparticle& superparticle,
                                            const Tendencies& tendencies) {
    state.change_layer(superparticle.z, state.grid, {0, 0, -tendencies.dqc, 0});
}

Tendencies ColumnModel::calc_tendencies(const Superparticle& superparticle,
                                        const double S, const double T,
                                        const double E, const double dt) {
    auto tendencies = condensation(superparticle.qc, superparticle.N,
                                   superparticle.r_dry, S, T, E, dt);
    return tendencies;
}

void ColumnModel::nucleation(Superparticle& s) {
    auto r = radius(s.qc, s.N, s.r_dry);
    if (std::abs((r - s.r_dry) / s.r_dry) < 1.e-5) {
        s.is_nucleated = false;
    }
}

bool ColumnModel::is_running() {
    runs++;
    state.t = runs * dt;
    if (state.t < t_max) {
        // if (runs < 12) {
        return true;
    } else {
        return false;
    }
}
