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
    virtual void setupdraft(State& state, double t){}
    virtual void keepcloudbase(State& state){};
};

inline void setupdraft(State& state, double t, double lifetime){
    for (unsigned int i = 0; i < state.levels.size(); ++i){
        state.levels[i].w = state.w_init * std::sin(2*PI*t/lifetime);
    }
}

inline void keepcloudbase(State& state, int n){
    int cloudbase_index = std::floor(state.cloud_base / state.grid.length) - 1;
    for(int i=0; i<n;++i){
        state.layers[cloudbase_index].qv = saturation_vapor(state.layers[cloudbase_index].T, 
                                                         state.layers[cloudbase_index].p);
        ++cloudbase_index;
    }
}

class AdvectFirstOrder : public Advect {
   public:
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

class AdvectAndSetFirstOrder: public AdvectFirstOrder{
    public:
        AdvectAndSetFirstOrder(double lifetime): lifetime(lifetime){}
        void setupdraft(State& state, double t){
            ::setupdraft(state, t, lifetime);
        }
    private:
        double lifetime;
};

class AdvectFirstOrderUpdraft : public AdvectAndSetFirstOrder {
    public:
        using AdvectAndSetFirstOrder::AdvectAndSetFirstOrder;
        void advect(State& state, const double& dt) override {
            auto cloud_bottom_lay = state.layers.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            auto cloud_bottom_lev = state.levels.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            first_order_upwind(member_iterator(cloud_bottom_lay, &Layer::qv),
                               member_iterator(state.layers.end(), &Layer::qv),
                               member_iterator(cloud_bottom_lev, &Level::w),
                               state.grid.length, dt);
        }
    private:
        double lifetime;
};

class AdvectSecondOrderUpdraft : public AdvectAndSetFirstOrder {
    public:
        using AdvectAndSetFirstOrder::AdvectAndSetFirstOrder;
        void advect(State& state, const double& dt) override {
            auto cloud_bottom_lay = state.layers.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            auto cloud_bottom_lev = state.levels.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            second_order_upwind(member_iterator(cloud_bottom_lay, &Layer::qv),
                               member_iterator(state.layers.end(), &Layer::qv),
                               member_iterator(cloud_bottom_lev, &Level::w),
                               state.grid.length, dt);
        }
        void keepcloudbase(State& state){
            ::keepcloudbase(state, 2);
        }
    private:
        double lifetime;
};

class AdvectSecondFirstOrderUpdraft : public AdvectAndSetFirstOrder {
    public:
        using AdvectAndSetFirstOrder::AdvectAndSetFirstOrder;
        void advect(State& state, const double& dt) override {
            auto cloud_bottom_lay = state.layers.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            auto cloud_bottom_lev = state.levels.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            second_first_order_upwind(member_iterator(cloud_bottom_lay, &Layer::qv),
                               member_iterator(state.layers.end(), &Layer::qv),
                               member_iterator(cloud_bottom_lev, &Level::w),
                               state.grid.length, dt);
        }
        void keepcloudbase(State& state){
            ::keepcloudbase(state, 2);
        }
    private:
        double lifetime;
};

class AdvectThirdOrderUpdraft : public AdvectAndSetFirstOrder {
    public:
        using AdvectAndSetFirstOrder::AdvectAndSetFirstOrder;
        void advect(State& state, const double& dt) override {
            auto cloud_bottom_lay = state.layers.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            auto cloud_bottom_lev = state.levels.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            third_order_upwind(member_iterator(cloud_bottom_lay, &Layer::qv),
                               member_iterator(state.layers.end(), &Layer::qv),
                               member_iterator(cloud_bottom_lev, &Level::w),
                               state.grid.length, dt);
        }
        void keepcloudbase(State& state){
            ::keepcloudbase(state, 3);
        }
    private:
        double lifetime;
};

class AdvectSixthOrderWickerSkamarock: public AdvectAndSetFirstOrder {
    public:
        using AdvectAndSetFirstOrder::AdvectAndSetFirstOrder;
        void advect(State& state, const double& dt) override {
            auto cloud_bottom_lay = state.layers.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            auto cloud_bottom_lev = state.levels.begin() +
                                    std::floor(state.cloud_base / state.grid.length) - 1;
            sixth_order_wickerskamarock(member_iterator(cloud_bottom_lay, &Layer::qv),
                               member_iterator(state.layers.end(), &Layer::qv),
                               member_iterator(cloud_bottom_lev, &Level::w),
                               state.grid.length, dt);
        }
        void keepcloudbase(State& state){
            ::keepcloudbase(state, 4);
        }
    private:
        double lifetime;
};


inline std::unique_ptr<Advect> mkAFO(){
    return std::make_unique<AdvectFirstOrder>();
}
inline std::unique_ptr<Advect> mkASFO(double lifetime){
    return std::make_unique<AdvectAndSetFirstOrder>(lifetime);
}
inline std::unique_ptr<Advect> mkAFOU(double lifetime){
    return std::make_unique<AdvectFirstOrderUpdraft>(lifetime);
}
inline std::unique_ptr<Advect> mkASOU(double lifetime){
    return std::make_unique<AdvectSecondOrderUpdraft>(lifetime);
}
inline std::unique_ptr<Advect> mkASFOU(double lifetime){
    return std::make_unique<AdvectSecondFirstOrderUpdraft>(lifetime);
}
inline std::unique_ptr<Advect> mkATOU(double lifetime){
    return std::make_unique<AdvectThirdOrderUpdraft>(lifetime);
}
inline std::unique_ptr<Advect> mkASOWK(double lifetime){
    return std::make_unique<AdvectSixthOrderWickerSkamarock>(lifetime);
}
