#include "collision.h"
#include "efficiencies.h"
#include "gtest/gtest.h"
#include <vector>
#include "superparticle.h"

TEST(hall_collision_kernal, test_value){
    double r = 1.e-6;
    double R = 1.e-6;
    double dfs = fall_speed(R)-fall_speed(r);
    HallCollisionKernal hk({});
    EXPECT_DOUBLE_EQ(hk(R, r, dfs), 0);
}

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
