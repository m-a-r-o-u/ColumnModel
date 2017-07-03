#include "gtest/gtest.h"
#include "state.h"
#include <algorithm>

TEST(state, swap_is_deep){
    State s0{0,{{0},1},{{0}, 1},{{0},1},{{0},1},{{0}, 1}};
    State s1{0,{{1},1},{{1}, 1},{{1},1},{{1},1},{{1}, 1}};
    std::swap(s0, s1);
    EXPECT_TRUE(s0.w(0.5, 1) == 1);
    EXPECT_TRUE(s1.w(0.5, 1) == 0);
}
