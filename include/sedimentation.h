#pragma once
#include "thermodynamic.h"
#include "twomey_utils.h"
#include "interpolate.h"
#include <vector>
#include <memory>

class Sedimentation{
    public:
    virtual double fall_speed(double r) const = 0;
 };

class FallSpeedLU: public Sedimentation{
    public:
    FallSpeedLU(){
        r = {0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4,
             2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8};
        std::transform(r.begin(), r.end(), r.begin(), [](double x){return x/2000.;});
        v = {0.27, 0.72, 1.17, 1.62, 2.06, 2.47, 2.87, 3.27, 3.67, 4.03, 4.64, 5.17, 5.65, 6.09, 6.49, 6.90, 7.27, 
             7.57, 7.82, 8.06, 8.26, 8.44, 8.60, 8.72, 8.83, 8.92, 8.98, 9.03, 9.07, 9.09, 9.12, 9.14, 9.16, 9.17};
    }
    double fall_speed(double rin) const override {
        double rout;
        if(rin<r[0]){
            double k1 = 1.19e8;
            rout = k1 * rin * rin;
        }
        else{
            int idx = left_index_min_zero_max_smallerlast(r, rin);
            rout = linear_interpolate(r[idx], v[idx], r[idx+1], v[idx+1], rin);
        }
        return std::min(rout, 10.);
    }
    private:
    std::vector<double> r;
    std::vector<double> v;
};

class NoFallSpeed: public Sedimentation {
    public:
    double fall_speed(double r) const override {
        return 0.;
    }
};

inline std::unique_ptr<FallSpeedLU> mkFSLU(){
    return std::make_unique<FallSpeedLU>();
}

inline std::unique_ptr<NoFallSpeed> mkNFS(){
    return std::make_unique<NoFallSpeed>();
}
