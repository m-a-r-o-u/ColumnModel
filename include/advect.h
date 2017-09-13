#pragma once
#include <cmath>
#include <memory>
#include <iostream>
#include "member_iterator.h"
#include "state.h"
#include "logger.h"

class Advect {
   public:
    virtual void advect(State& state, const double& dt) = 0;
    virtual void init(Logger& logger){}
};

class AdvectFirstOrder : public Advect {
   public:
    AdvectFirstOrder() {}
    void advect(State& state, const double& dt) override {
        auto cloud_bottom_lay = state.layers.begin() +
                                std::floor(state.cloud_base / state.grid.length) - 1;
        auto cloud_bottom_lev = state.levels.begin() +
                                std::floor(state.cloud_base / state.grid.length) - 1;
        advect_first_order(member_iterator(cloud_bottom_lay, &Layer::qv),
                           member_iterator(state.layers.end(), &Layer::qv),
                           member_iterator(cloud_bottom_lev, &Level::w),
                           state.grid.length, dt);
    }
};

 inline std::unique_ptr<Advect> mkFirstOrder(){
     return std::make_unique<AdvectFirstOrder>();
 }

