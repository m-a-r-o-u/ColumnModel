#pragma once
#include <memory>
#include "state.h"
#include "tendencies.h"
#include "superparticle.h"
#include "superparticle_source.h"
#include "logger.h"
#include "superparticle_population.h"
#include "grid.h"

class ColumnModel {
   public:
    typedef std::back_insert_iterator<std::vector<Superparticle>> OIt;
    ColumnModel(const State& initial_state,
                std::shared_ptr<SuperParticleSource<OIt>> source, double t_max,
                double dt, double superparticle_dz, Grid grid)
        : state(initial_state),
          source(source),
          t_max(t_max),
          dt(dt),
          superparticles{{}, superparticle_dz},
          grid(grid){
          };
    void run(std::shared_ptr<Logger> logger);

   private:
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
                               const double S, const double T,
                               const double E, const double dt);
    std::shared_ptr<SuperParticleSource<OIt>> source;
    const Grid grid;
    State state;
    Superparticles superparticles;
    const double dt, t_max;
    int runs = 0;
};
