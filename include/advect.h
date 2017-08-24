#pragma once
#include <cmath>
#include <memory>
#include "member_iterator.h"
#include "state.h"

class Advect {
   public:
    virtual void advect(State& state, const double& dt) = 0;
};

class AdvectFirstOrder : public Advect {
   public:
    AdvectFirstOrder(const double& cloud_base): cloud_base(cloud_base){}
    void advect(State& state, const double& dt) override {
        auto cloud_bottom_lay = state.layers.begin() +
                                std::floor(cloud_base / state.grid.length) - 1;
        auto cloud_bottom_lev = state.levels.begin() +
                                std::floor(cloud_base / state.grid.length) - 1;
        advect_first_order(member_iterator(cloud_bottom_lay, &Layer::qv),
                           member_iterator(state.layers.end(), &Layer::qv),
                           member_iterator(cloud_bottom_lev, &Level::w),
                           state.grid.length, dt);
    }

   private:
    double cloud_base;
};

 inline std::unique_ptr<Advect> mkFirstOrder(const double& cloud_base){
     return std::make_unique<AdvectFirstOrder>(cloud_base);
 }

