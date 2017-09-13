#pragma once
#include <memory>
#include "advect.h"
#include "grid.h"
#include "logger.h"
#include "radiationsolver.h"
#include "saturation_fluctuations.h"
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
                std::unique_ptr<FluctuationSolver> fluctuations)
        : source(source),
          state(initial_state),
          superparticles{},
          dt(dt),
          t_max(t_max),
          radiation_solver(radiation_solver),
          grid(std::move(grid)),
          advection_solver(std::move(advection_solver)),
          fluctuations(std::move(fluctuations)){};
    void run(std::shared_ptr<Logger> logger);

   private:
    void log_every_seconds(std::shared_ptr<Logger> logger, double dt_out);
    void step();
    bool is_running();
    void apply_tendencies_to_superparticle(Superparticle& superparticle,
                                           Tendencies& tendencies);
    void apply_tendencies_to_state(Superparticle& superparticle,
                                   const Tendencies& tendencies);
    void insert_superparticles();
    void nucleation(Superparticle& superparticle);
    Tendencies calc_tendencies(const Superparticle& superparticle,
                               const double S, const double T, const double E,
                               const double dt);
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
};
