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
    std::vector<double> cooling(state.grid.n_lay, - 2.3e-5 / dt);
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

void check_superparticles(std::vector<Superparticle>& sp, const Grid& grid) {
    for (auto s : sp) {
        if (s.N > sp[0].N || s.N < sp[0].N) {
            std::cout << "somethin with the multiplicity" << std::endl;
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
        log_every_seconds(logger, 30.);
    }
}

void ColumnModel::step() {
    source->generateParticles(std::back_inserter(superparticles), state, dt,
                                  superparticles);
    advection_solver->advect(state, dt);
    fluctuations->refresh(superparticles);

    State old_state(state);

    if (true) {
        check_state(state);
        check_superparticles(superparticles, state.grid);
    }

    for (auto& sp : superparticles) {
        Layer lay = old_state.layer_at(sp.z);
        Level lvl = old_state.upper_level_at(sp.z);

        double S = super_saturation(lay.T, lay.p, lay.qv) + fluctuations->getFluctuation(sp, dt);

        if (sp.is_nucleated) {
            auto tendencies = calc_tendencies(sp, S, lay.T, lay.E, dt);
            apply_tendencies_to_state(sp, tendencies);
            apply_tendencies_to_superparticle(sp, tendencies, lvl);
        }
    }
    removeUnnucleated(superparticles);
    radiation_solver.calculate_radiation(state, superparticles);

    cooling_the_column(state, dt);
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
    Superparticle& sp, Tendencies& tendencies,
    const Level& lvl) {
    sp.z += lvl.w * dt - dt * 0.1;//fall_speed(radius(sp.qc, sp.N, sp.r_dry));
    sp.qc += tendencies.dqc;
    nucleation(sp);
}

void ColumnModel::apply_tendencies_to_state(const Superparticle& sp,
                                            const Tendencies& tendencies) {
    state.change_layer(sp.z, {0, 0, -tendencies.dqc, 0});
}

Tendencies ColumnModel::calc_tendencies(const Superparticle& superparticle,
                                        const double S, const double T,
                                        const double E, const double dt) {
    auto tendencies = condensation(superparticle.qc, superparticle.N,
                                   superparticle.r_dry, S, T, E, dt);
    return tendencies;
}

void ColumnModel::nucleation(Superparticle& s) {
    if (s.qc <= 0) {
        s.is_nucleated = false;
    }
    if (s.z <= 0){
        s.is_nucleated = false;
    }
//    if (s.radius() >= 20.e-6){
//        s.is_nucleated = false;
//    }
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
