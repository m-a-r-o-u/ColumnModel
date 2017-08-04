#pragma once
#include <cmath>
#include "linearfield.h"
#include "layer_quantities.h"
#include "level_quantities.h"
#include "grid.h"

struct State {
    double t;
    std::vector<Layer> layers;
    std::vector<Level> levels;

    inline Layer& layer_at(double z, Grid grid) {
        int index = std::floor(z / grid.length);
        return layers[index];
    }
    inline void change_layer(double z, Grid grid, const Layer& tendencies){
        int index = std::floor(z / grid.length);
        layers[index] += tendencies;
    }
    inline const Level& lower_level_at(double z, Grid grid) {
        int index = std::floor(z / grid.length);
        return levels[index];
    }
    inline const Level& upper_level_at(double z, Grid grid) {
        int index = std::ceil(z / grid.length);
        return levels[index];
    }
};
