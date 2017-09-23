#include <vector>
#include <fstream>
#include <iostream>
#include "grid.h"
#include "gtest/gtest.h"
#include "saturation_fluctuations.h"
#include "superparticle.h"

TEST(saturation_fluctuations, test_tke){
    double EPSILON = 1.e-2;
    double l = 10.;
    double e = 10.e-4;
    ASSERT_TRUE(std::abs(turbulent_kinetic_energy(l, e) - 5.2e-2) < EPSILON);
}

TEST(saturation_fluctuations, test_integral_time_scale){
    double EPSILON = 0.2;
    double l = 10.;
    double e = 10.e-4;
    double tke = turbulent_kinetic_energy(l, e);
    ASSERT_TRUE(std::abs(integral_timescale(l, tke) - 29.) < EPSILON);
}

TEST(saturation_fluctuations, test_ornstein_uhlenbeck_process){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    double w = 0.;
    double w_max(w);
    double dt = 0.02;
    double w_std = 0.5440625136744149;
    double tau = 49.80365826707145;
    int i_end = 20 * 60 / dt;

    for (int t = 0; t < i_end;++t){
        w = ornstein_uhlenbeck_process(gen, w, dt, tau, w_std);
        if (std::abs(w) > std::abs(w_max)){w_max = w;}
        ASSERT_LT(std::abs( w_max ), 5.);
    }
}

TEST(saturation_fluctuations, test_saturation_fluctuations){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    Grid  grid{300., 100.};

    std::vector<Superparticle> sp;
    sp.push_back({0.00001, 50, 1.e-6, 100000000, true});

    auto fsolver = mkFS(gen, "markov", 50.e-4, grid);

    double dt = 0.1;
    double tmax = 20.*60. * 10.;

    std::ofstream mfile;
    mfile.open("./test/test_saturation_fluctuations.txt");
    for (int t = 0; t< tmax;++t){
        fsolver->refresh(sp);
        fsolver->getFluctuation(sp[0], dt);
        mfile << sp[0].S_prime << " " << sp[0].w_prime << "\n";
    }
    mfile.close();
}
