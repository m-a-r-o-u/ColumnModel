#include <numeric>
#include <vector>
#include "gtest/gtest.h" 

TEST(accumulate, test_ints){
     std::vector<int> x{1,2,3};
     int sum = std::accumulate(x.begin(), x.end(), 0);
     ASSERT_TRUE(sum == 6);
}
TEST(accumulate, test_doubles){
     std::vector<double> x{1.,2.,3.};
     double sum = std::accumulate(x.begin(), x.end(), 0.);
     ASSERT_TRUE(sum == 6.);
}
