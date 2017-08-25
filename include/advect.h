#pragma once
#include <cmath>
#include <memory>
#include "member_iterator.h"
#include "state.h"

class Advect {
   public:
    virtual void advect(State& state, const double& dt) = 0;
    void init(Logger& logger){}
};

class AdvectFirstOrder : public Advect {
   public:
    AdvectFirstOrder(const double& cloud_base, const double& w_init): cloud_base(cloud_base), w_init(w_init){}
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
    void init(Logger& logger){
        logger.setAttr("w_init", w_init);
    }

   private:
    double cloud_base;
    double w_init;
};

 inline std::unique_ptr<Advect> mkFirstOrder(const double& cloud_base, const double& w_init){
     return std::make_unique<AdvectFirstOrder>(cloud_base, w_init);
 }

