#include <grid.h>
#include <vector>
#include "analize_sp.h"
#include "gtest/gtest.h"
#include "superparticle.h"
#include "thermodynamic.h"

TEST(count_nucleated, all_are_nucleated) {
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), true},
                                 {0.00001, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<int> res = count_nucleated(v, grid);

    EXPECT_TRUE(res == std::vector<int>({0, 2, 1}));
}

TEST(count_nucleated, not_all_are_nucleated) {
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<int> res = count_nucleated(v, grid);

    EXPECT_TRUE(res == std::vector<int>({0, 1, 1}));
}

TEST(calculate_qc_profile, not_all_are_nucleated) {
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_qc_profile(v, grid);

    EXPECT_TRUE(res == std::vector<double>({0, 1.e-5, 1.e-5}));
}

TEST(calculate_effective_radius_profile, check_values) {
    double EPSILON = 1.e-5;
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_effective_radius_profile(v, grid);
    std::vector<double> cp(
        {0, radius(1.e-5, 1e8, 1.e-6), radius(1.e-5, 1e8, 1.e-6)});
    for (unsigned int i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(std::abs(res[i] - cp[i]) < EPSILON);
    }
}

TEST(calculate_maximal_radius_profile, check_values) {
    double EPSILON = 1.e-5;
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), true},
                                 {0.00002, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_maximal_radius_profile(v, grid);
    std::vector<double> cp(
        {0, radius(1.e-5, 1e8, 1.e-6), radius(2.e-5, 1e8, 1.e-6)});
    for (unsigned int i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(std::abs(res[i] - cp[i]) < EPSILON);
    }
}

TEST(calculate_minimal_radius_profile, check_values) {
    double EPSILON = 1.e-5;
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), true},
                                 {0.00002, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_minimal_radius_profile(v, grid);
    std::vector<double> cp(
        {0, radius(1.e-5, 1e8, 1.e-6), radius(1.e-5, 1e8, 1.e-6)});
    for (unsigned int i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(std::abs(res[i] - cp[i]) < EPSILON);
    }
}

TEST(calculate_mean_radius_profile, check_values) {
    double EPSILON = 1.e-5;
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), true},
                                 {0.00001, 2, 1.e-6, int(1e8), true},
                                 {0.00001, 2, 1.e-6, int(1e8), true},
                                 {0.00001, 2, 1.e-6, int(1e8), true}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_mean_radius_profile(v, grid);
    std::vector<double> cp(
        {0, radius(1.e-5, 1e8, 1.e-6), radius(1.e-5, 1e8, 1.e-6)});
    for (unsigned int i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(std::abs(res[i] - cp[i]) < EPSILON);
    }
}

TEST(calculate_stddev_radius_profile, no_nucleated) {
    double EPSILON = 1.e-5;
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), false},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), false}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_stddev_radius_profile(v, grid);
    std::vector<double> cp(
        {0, 0, 0});
    for (unsigned int i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(std::abs(res[i] - cp[i]) < EPSILON);
    }
}

TEST(calculate_stddev_radius_profile, one_nucleated) {
    double EPSILON = 1.e-5;
    std::vector<Superparticle> v{{0.00001, 1, 1.e-6, int(1e8), true},
                                 {0.00001, 1.4, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), false},
                                 {0.00001, 2, 1.e-6, int(1e8), false}};
    Grid grid{3., 1.};
    std::vector<double> res = calculate_stddev_radius_profile(v, grid);
    std::vector<double> cp(
        {0, 0, 0});
    for (unsigned int i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(std::abs(res[i] - cp[i]) < EPSILON);
    }
}
