#include "collision.h"
#include "gtest/gtest.h"
#include <vector>
#include "superparticle.h"
#include "efficiencies.h"

//TEST(collision_efficiencies, test_low_left){
//    double EPSI = 1.e-10;
//    double R = 250.;
//    double rR = 0.075;
//    double test = ((0.87 + 0.96) + (1.00 + 0.97)) / 4.;
//    EXPECT_TRUE(std::abs(test - collision_efficiency(R, rR)) < EPSI);
//}
//TEST(collision_efficiencies, test_low_right){
//    double EPSI = 1.e-10;
//    double R = 250.;
//    double rR = 0.975;
//    double test = 1.00;
//    EXPECT_TRUE(std::abs(test - collision_efficiency(R, rR)) < EPSI);
//}
//
//TEST(collision_efficiencies, test_up_right){
//    double EPSI = 1.e-10;
//    double R = 15.;
//    double rR = 0.975;
//    double test = (0.029+0.027 + 0.033 + 0.027) / 4.;
//    EXPECT_TRUE(std::abs(test - collision_efficiency(R, rR)) < EPSI);
//}
//TEST(collision_efficiencies, test_up_left){
//    double EPSI = 1.e-10;
//    double R = 15.;
//    double rR = 0.075;
//    double test = 0.0001;
//    EXPECT_TRUE(std::abs(test - collision_efficiency(R, rR)) < EPSI);
//}
//TEST(hall_collision_kernal, test_value){
//    double EPSI = 1.e-10;
//    double r = 1.e-6;
//    double R = 1.e-6;
//    EXPECT_TRUE(std::abs(0 - hall_collision_kernal(R, r)) < EPSI);
//}
//TEST(collide, test_two_values){
//    double EPSI = 1.e-10;
//    std::vector<Superparticle> sps{{0.001, 0, 0, 100000000, true},
//                                   {0.003, 0, 0, 100000002, true}
//                                   };
//    double dt = 0.1;
//    auto out = collide(sps, dt);
//    EXPECT_TRUE(std::abs(out[0].dqc + out[1].dqc) < EPSI);
//    EXPECT_TRUE(out[0].dN <= 0);
//    EXPECT_TRUE(out[1].dN == 0);
//}
//TEST(collide, test_tree_values){
//    double EPSI = 1.e-10;
//    std::vector<Superparticle> sps{{0.001, 0, 0, 100000001, true},
//                                   {0.002, 0, 0, 100000002, true},
//                                   {0.003, 0, 0, 100000003, true}
//                                   };
//    double dt = 0.1;
//    auto out = collide(sps, dt);
//    EXPECT_TRUE(std::abs(out[0].dqc + out[1].dqc + out[2].dqc) < EPSI);
//    EXPECT_TRUE(out[0].dN <= 0);
//    EXPECT_TRUE(out[2].dN == 0);
//}
TEST(sortsps, test_nvalues){
    Grid grid(1000., 250.);
    std::vector<Superparticle> sps{
                                   {0.005, 125, 0, 100000000, true},
                                   {0.001, 375, 0, 100000000, true},
                                   {0.005, 625, 0, 100000000, true},
                                   {0.005, 875, 0, 100000000, true},
                                   {0.001, 625, 0, 100000000, true},
                                   {0.001, 125, 0, 100000000, true},
                                   {0.005, 375, 0, 100000000, true},
                                   {0.001, 875, 0, 100000000, true},
                                   };
    auto out = sortsps(sps, grid);
    EXPECT_TRUE(out[0][0].i == 5); 
    EXPECT_TRUE(out[0][1].i == 0); 
    EXPECT_TRUE(out[1][0].i == 1); 
    EXPECT_TRUE(out[1][1].i == 6); 
    EXPECT_TRUE(out[2][0].i == 4); 
    EXPECT_TRUE(out[2][1].i == 2); 
    EXPECT_TRUE(out[3][0].i == 7); 
    EXPECT_TRUE(out[3][1].i == 3); 
}
//
//TEST(sortsps, test_particle_in_ground){
//    Grid grid(1000., 250.);
//    std::vector<Superparticle> sps{
//                                  {-0.005, 1, 0, 100000000, false},
//                                  };
//    auto out = sortsps(sps, grid);
//    EXPECT_TRUE(out[0][0].i == 5); 
//}
