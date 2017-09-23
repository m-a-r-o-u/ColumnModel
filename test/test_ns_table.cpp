#include "ns_table.h"
#include "gtest/gtest.h" 
#include <iostream>
#include <algorithm>

TEST(load_data, test_is_sorted){
    std::vector<double> n;
    std::vector<double> s;
    load_data(n, s);
    EXPECT_TRUE(std::is_sorted(n.begin(), n.end()) == true);
    EXPECT_TRUE(std::is_sorted(s.begin(), s.end()) == true);
}

TEST(load_data, test_starts_zero){
    std::vector<double> n;
    std::vector<double> s;
    load_data(n, s);
    EXPECT_TRUE(n[0]==0);
    EXPECT_TRUE(s[0]==0);
}

TEST(arange, test_ints){
    int start = 0;
    int stop = 2;
    int step = 1;
    auto out = arange(start, stop, step);
    EXPECT_TRUE(out[0]==0);
    EXPECT_TRUE(out[1]==1);
    EXPECT_TRUE(out.size() == 2);
}

TEST(arange, test_doubles){
    double start = 0;
    double stop = 2;
    double step = 1;
    auto out = arange(start, stop, step);
    EXPECT_TRUE(out[0]==0.);
    EXPECT_TRUE(out[1]==1.);
    EXPECT_TRUE(out.size() == 2);
}
TEST(arange, test_mix_ints_doubles){
    double start = 0;
    double stop = 2;
    int step = 1;
    auto out = arange(start, stop, step);
    EXPECT_TRUE(out[0]==0.);
    EXPECT_TRUE(out[1]==1.);
    EXPECT_TRUE(out.size() == 2);
}
TEST(arange, test_stop_lower_start){
    double start = 0;
    double stop = -1;
    int step = 1;
    auto out = arange(start, stop, step);
    EXPECT_TRUE(out.size()==0.);
}

TEST(arange, test_negative_step){
    double start = 0;
    double stop = 2;
    int step = -1;
    auto out = arange(start, stop, step);
    EXPECT_TRUE(out.size()==0.);
}

TEST(linear_interpolate, test_value){
    double x1 = 1;
    double y1 = 1;
    double x2 = 3;
    double y2 = 3;
    double x = 2;
    auto out = linear_interpolate(x1,y1,x2,y2,x);
    EXPECT_TRUE(out == 2);
}

TEST(linear_interpolate, test_value_on_value){
    double x1 = 1;
    double y1 = 1;
    double x2 = 3;
    double y2 = 3;
    double x = 1;
    auto out = linear_interpolate(x1,y1,x2,y2,x);
    EXPECT_TRUE(out == 1);
}
TEST(linear_interpolate, test_value_on_value_2){
    double x1 = 1;
    double y1 = 1;
    double x2 = 3;
    double y2 = 3;
    double x = 3;
    auto out = linear_interpolate(x1,y1,x2,y2,x);
    EXPECT_TRUE(out == 3);
}

TEST(linear_interpolate, test_smaller_then_min){
    double x1 = 1;
    double y1 = 1;
    double x2 = 3;
    double y2 = 3;
    double x = 0;
    auto out = linear_interpolate(x1,y1,x2,y2,x);
    EXPECT_TRUE(out == 0);
}

TEST(linear_interpolate, test_larger_then_max){
    double x1 = 1;
    double y1 = 1;
    double x2 = 3;
    double y2 = 3;
    double x = 4;
    auto out = linear_interpolate(x1,y1,x2,y2,x);
    EXPECT_TRUE(out == 4);
}
