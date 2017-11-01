#include "gtest/gtest.h"
#include "projection_iterator.h"
#include <vector>

TEST(projection_iterator, test_value){
    std::vector<std::vector<int>> foo{{1,2},{3,4}};
    auto it = projection_iterator(foo.begin(), [](auto x){ return x[0];});
    EXPECT_EQ(*it, 1);
}

TEST(projection_iterator, test_write_value){
    std::vector<std::vector<int>> foo{{1,2},{3,4}};
    auto it = projection_iterator(foo.begin(), [](auto& x) -> int& { return x[0];});
    *it = 5;
    EXPECT_EQ(foo[0][0], 5);
}

