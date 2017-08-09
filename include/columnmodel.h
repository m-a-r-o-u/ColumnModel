#pragma once
#include <memory>
#include "grid.h"
#include "logger.h"
#include "radiationsolver.h"
#include "state.h"
#include "superparticle.h"
#include "superparticle_source.h"
#include "tendencies.h"

class ColumnModel {
   public:
    typedef std::back_insert_iterator<std::vector<Superparticle>> OIt;
    ColumnModel(const State& initial_state,
                std::shared_ptr<SuperParticleSource<OIt>> source, double t_max,
                double dt, Grid grid, RadiationSolver radiation_solver, bool sw,
                bool lw)
        : source(source),
          grid(grid),
          state(initial_state),
          superparticles{},
          dt(dt),
          t_max(t_max),
          radiation_solver(radiation_solver),
          lw(lw),
          sw(sw)
          {};
    void run(std::shared_ptr<Logger> logger);

   private:
    void log_every_seconds(std::shared_ptr<Logger> logger, double dt_out);
    void step();
    bool is_running();
    void apply_tendencies_to_superparticle(Superparticle& superparticle,
                                           Tendencies& tendencies,
                                           const double w);
    void apply_tendencies_to_state(Superparticle& superparticle,
                                   const Tendencies& tendencies);
    void insert_superparticles();
    void nucleation(Superparticle& superparticle);
    Tendencies calc_tendencies(const Superparticle& superparticle,
                               const double S, const double T, const double E,
                               const double dt);
    std::shared_ptr<SuperParticleSource<OIt>> source;
    const Grid grid;
    State state;
    std::vector<Superparticle> superparticles;
    const double dt;
    const double t_max;
    int runs = 0;
    RadiationSolver radiation_solver;
    bool lw;
    bool sw;
};
