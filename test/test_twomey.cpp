#include "twomey.h"
#include "twomey_utils.h"
#include "gtest/gtest.h" 
#include <iostream>

TEST(test_lower_bound, test_xsmaller){
    std::vector<double> r{1,2,3,4};
    double x = 0;
    auto itr = std::lower_bound(r.begin(), r.end(), x);
    ASSERT_TRUE(*itr == 1);
}
TEST(test_lower_bound, test_xlarger){
    std::vector<double> r{1,2,3,4};
    double x = 5;
    auto itr = std::lower_bound(r.begin(), r.end(), x);
    ASSERT_TRUE(itr == r.end());
}
TEST(test_lower_bound, test_x_in_r){
    std::vector<double> r{1,2,3,4};
    double x = 2.5;
    auto itr = std::lower_bound(r.begin(), r.end(), x);
    ASSERT_TRUE(*itr == 3);
}
TEST(test_lower_bound, test_x_element_of_r){
    std::vector<double> r{1,2,3,4};
    double x = 2;
    auto itr = std::lower_bound(r.begin(), r.end(), x);
    ASSERT_TRUE(*itr == 2);
}

TEST(test_lower_bound_index, test_x_smaller){
    std::vector<double> r {1.,2.,3.};
    double x = 0;
    auto i = lower_bound_index(r.begin(), r.end(), x);
    ASSERT_TRUE(i == 0);
}
TEST(test_lower_bound_index, test_x_larger){
    std::vector<double> r {1.,2.,3.};
    double x = 4;
    auto i = lower_bound_index(r.begin(), r.end(), x);
    ASSERT_TRUE(i == 3);
}
TEST(test_lower_bound_index, test_x_in_range){
    std::vector<double> r {1.,2.,3.};
    double x = 1.5;
    auto i = lower_bound_index(r.begin(), r.end(), x);
    ASSERT_TRUE(i == 1);
}
TEST(test_lower_bound_index, test_x_element_of_r){
    std::vector<double> r {1.,2.,3.};
    double x = 2;
    auto i = lower_bound_index(r.begin(), r.end(), x);
    ASSERT_TRUE(i == 1);
}
TEST(test_upper_bound_index, test_x_element_of_r){
    std::vector<double> r {1.,2.,3.};
    double x = 2;
    auto i = upper_bound_index(r.begin(), r.end(), x);
    ASSERT_TRUE(i == 2);
}
TEST(indexes, test_a_range_of_elements){
    std::vector<double> r {1.,2.,3.};
    std::vector<double> x {-1, 1, 1.5, 2, 4};
    auto i = indexes(r, x);
    ASSERT_TRUE(i[0] == 0);
    ASSERT_TRUE(i[1] == 1);
    ASSERT_TRUE(i[2] == 1);
    ASSERT_TRUE(i[3] == 2);
    ASSERT_TRUE(i[4] == 3);
}

TEST(left_index, test_try_value){
    std::vector<double> r {1.,2.,3.};
    double x = 1.2;
    auto i = left_index(r, x);
    ASSERT_TRUE(i == 0);
}

TEST(left_index, test_try_value_2){
    std::vector<double> r {1.,2.,3.};
    double x = 2.;
    auto i = left_index(r, x);
    ASSERT_TRUE(i == 0);
}
TEST(left_index, test_larger_then_max){
    std::vector<double> r {1.,2.,3.};
    double x = 100.;
    auto i = left_index(r, x);
    ASSERT_TRUE(i == 2);
}

TEST(left_index, test_smaller_min){
    std::vector<double> r {1.,2.,3.};
    double x = 1.;
    auto i = left_index(r, x);
    ASSERT_TRUE(i == -1);
}

TEST(right_index, test_try_value){
    std::vector<double> r {1.,2.,3.};
    double x = 1.2;
    auto i = right_index(r, x);
    ASSERT_TRUE(i == 1);
}

TEST(right_index, test_try_value_2){
    std::vector<double> r {1.,2.,3.};
    double x = 2.;
    auto i = right_index(r, x);
    ASSERT_TRUE(i == 1);
}
TEST(right_index, test_larger_then_max){
    std::vector<double> r {1.,2.,3.};
    double x = 100.;
    auto i = right_index(r, x);
    ASSERT_TRUE(i == 3);
}

TEST(right_index, test_smaller_min){
    std::vector<double> r {1.,2.,3.};
    double x = 1.;
    auto i = right_index(r, x);
    ASSERT_TRUE(i == 0);
}
TEST(left_index_min_zero, test_try_value){
    std::vector<double> r {1.,2.,3.};
    double x = 1.2;
    auto i = left_index_min_zero(r, x);
    ASSERT_TRUE(i == 0);
}

TEST(left_index_min_zero, test_try_value_2){
    std::vector<double> r {1.,2.,3.};
    double x = 2.;
    auto i = left_index_min_zero(r, x);
    ASSERT_TRUE(i == 0);
}
TEST(left_index_min_zero, test_larger_then_max){
    std::vector<double> r {1.,2.,3.};
    double x = 100.;
    auto i = left_index_min_zero(r, x);
    ASSERT_TRUE(i == 2);
}

TEST(left_index_min_zero, test_smaller_min){
    std::vector<double> r {1.,2.,3.};
    double x = 1.;
    auto i = left_index_min_zero(r, x);
    ASSERT_TRUE(i == 0);
}
TEST(left_index_min_zero_max_smallerlast, test_larger_max){
    std::vector<double> r {1.,2.,3.};
    double x = 4.;
    auto i = left_index_min_zero_max_smallerlast(r, x);
    ASSERT_TRUE(i == 1);
}
