#pragma once
#include "state.h"
#include "level_quantities.h"
#include "layer_quantities.h"

class DummyContinuousStateView {
   public:
    inline DummyContinuousStateView(State& state, const Grid& grid)
        : state(state), grid(grid) {}

    inline Layer layer_at(double z, double dz) const {
        return state.layers[z / grid.length];
    }
    inline Level level_at(double z, double dz) const {
        return state.levels[z / grid.length];
    }
    inline void layer_change(double z, double dz, const Layer& val) {
        state.layers[z / grid.length] += val;
    }
    inline void level_change(double z, double dz, const Level& val) {
        state.levels[z / grid.length] += val;
    }

   private:
    State& state;
    const Grid& grid;
};
