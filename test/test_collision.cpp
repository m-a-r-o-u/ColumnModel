#include "collision.h"
#include "sedimentation.h"
#include "efficiencies.h"
#include "gtest/gtest.h"
#include <vector>
#include "superparticle.h"

TEST(hall_collision_kernal, test_value){
    double r = 1.e-6;
    double R = 1.e-6;
    double dfs = fall_speed(R)-fall_speed(r);
    HallCollisionKernal<Efficiencies> hk({});
    EXPECT_DOUBLE_EQ(hk(R, r, dfs), 0);
}

inline std::ostream& operator<<(std::ostream& str, const SpMassTendencies& mt) {
    str << "MT(" << mt.dN << ", " << mt.dqc << ")";
    return str;
}

TEST(collide, test_two_values){
    std::vector<Superparticle> sps{{0.001, 0, 0, 100000000, true},
                                   {0.003, 0, 0, 100000002, true}
                                   };
    double dt = 0.1;
    FallSpeedLU sedi;
    BoxCollisions<HallCollisionKernal<Efficiencies>> bc(sedi,{{}});
    std::vector<SpMassTendencies> mt(sps.size());
    bc.collide(sps.begin(), sps.end(), mt.begin(), dt);
    EXPECT_NEAR(mt[0].dqc + mt[1].dqc, 0, 1e-15);
    EXPECT_LE(mt[0].dN, 0);
    EXPECT_EQ(mt[1].dN, 0);
}

TEST(collide, test_two_values_unitefficiencies){
    std::vector<Superparticle> sps{{0.001, 0, 0, 100000000, true},
                                   {0.003, 0, 0, 100000002, true}
                                   };
    double dt = 0.1;
    FallSpeedLU sedi;
    BoxCollisions<HallCollisionKernal<UnitEfficiencies>> bc(sedi,{{}});
    std::vector<SpMassTendencies> mt(sps.size());
    bc.collide(sps.begin(), sps.end(), mt.begin(), dt);
    EXPECT_NEAR(mt[0].dqc + mt[1].dqc, 0, 1e-15);
    EXPECT_LE(mt[0].dN, 0);
    EXPECT_EQ(mt[1].dN, 0);
}
TEST(collide, test_same_two_values){
    std::vector<Superparticle> sps{{0.001, 0, 0, 100000000, true},
                                   {0.001, 0, 0, 100000000, true}
                                   };
    double dt = 0.1;
    FallSpeedLU sedi;
    BoxCollisions<HallCollisionKernal<Efficiencies>> bc(sedi,{{}});
    std::vector<SpMassTendencies> mt(sps.size());
    bc.collide(sps.begin(), sps.end(), mt.begin(), dt);
    EXPECT_NEAR(mt[0].dqc + mt[1].dqc, 0, 1e-15);
    EXPECT_EQ(mt[0].dN, 0);
    EXPECT_EQ(mt[1].dN, 0);
}
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
