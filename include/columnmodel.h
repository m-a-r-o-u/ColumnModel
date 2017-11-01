#pragma once
#include <cstdlib>
#include <memory>
#include "advect.h"
#include "collision.h"
#include "grid.h"
#include "logger.h"
#include "radiationsolver.h"
#include "saturation_fluctuations.h"
#include "sedimentation.h"
#include "state.h"
#include "superparticle.h"
#include "superparticle_source.h"
#include "tendencies.h"

class ColumnModel {
   public:
    typedef std::back_insert_iterator<std::vector<Superparticle>> OIt;
    ColumnModel(const State& initial_state,
                std::shared_ptr<SuperParticleSource<OIt>> source, double t_max,
                double dt, RadiationSolver radiation_solver,
                std::unique_ptr<Grid> grid,
                std::unique_ptr<Advect> advection_solver,
                std::unique_ptr<FluctuationSolver> fluctuations,
                std::unique_ptr<Collisions> collisions,
                std::unique_ptr<Sedimentation> sedimentation)
        : source(source),
          state(initial_state),
          superparticles{},
          dt(dt),
          t_max(t_max),
          radiation_solver(radiation_solver),
          grid(std::move(grid)),
          advection_solver(std::move(advection_solver)),
          fluctuations(std::move(fluctuations)),
          collisions(std::move(collisions)),
          sedimentation(std::move(sedimentation)){};
    void run(std::shared_ptr<Logger> logger);

   private:
    void log_every_seconds(std::shared_ptr<Logger> logger, double dt_out);
    void step();
    bool is_running();
    void apply_tendencies_to_superparticle(Superparticle& superparticle,
                                           Tendencies& tendencies,
                                           const Level& lvl);
    void apply_tendencies_to_state(const Superparticle& superparticle,
                                   const Tendencies& tendencies);
    void apply_collision_tendencies(
        std::vector<Superparticle>& sps,
        const std::vector<SpMassTendencies>& tendencies);
    void insert_superparticles();
    Tendencies calc_tendencies(const Superparticle& superparticle,
                               const double S, const double T, const double E,
                               const double dt);

    void do_condensation(State& old_state);
    void do_collisions();
    std::shared_ptr<SuperParticleSource<OIt>> source;
    State state;
    std::vector<Superparticle> superparticles;
    const double dt;
    const double t_max;
    int runs = 0;
    RadiationSolver radiation_solver;
    std::unique_ptr<Grid> grid;
    std::unique_ptr<Advect> advection_solver;
    std::unique_ptr<FluctuationSolver> fluctuations;
    std::unique_ptr<Collisions> collisions;
    std::unique_ptr<Sedimentation> sedimentation;
};
