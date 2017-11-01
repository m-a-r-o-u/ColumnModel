#include "gtest/gtest.h"
#include <vector>
#include "superparticle.h"
#include "efficiencies.h"

TEST(collision_efficiencies, test_low_left){
    Efficiencies effi;
    double R = 250.;
    double rR = 0.075;
    double test = ((0.87 + 0.96) + (1.00 + 0.97)) / 4.;
    EXPECT_DOUBLE_EQ(effi.collision_efficiency(R, rR), test);
}
TEST(collision_efficiencies, test_low_right){
    Efficiencies effi;
    double R = 250.;
    double rR = 0.975;
    double test = 1.00;
    EXPECT_DOUBLE_EQ(effi.collision_efficiency(R, rR), test);
}

TEST(collision_efficiencies, test_up_right){
    Efficiencies effi;
    double R = 15.;
    double rR = 0.975;
    double test = (0.029+0.027 + 0.033 + 0.027) / 4.;
    EXPECT_DOUBLE_EQ(effi.collision_efficiency(R, rR), test);
}
TEST(collision_efficiencies, test_up_left){
    Efficiencies effi;
    double R = 15.;
    double rR = 0.075;
    double test = 0.0001;
    EXPECT_DOUBLE_EQ(effi.collision_efficiency(R, rR), test);
}
