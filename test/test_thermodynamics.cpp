#include "gtest/gtest.h"
#include "thermodynamic.h"

TEST(will_nucleate, s_smaller_zero) {
    EXPECT_FALSE(will_nucleate(1.e-6, -0.01, 280));
}

TEST(will_nucleate, s_larger_zero_but_smaller_s_critical) {
    EXPECT_FALSE(will_nucleate(1.e-8, 0.01, 273.15));
}

TEST(will_nucleate, s_larger_zero_and_larger_s_critical) {
    EXPECT_TRUE(will_nucleate(1.e-6, 0.01, 273.15));
}

TEST(critical_supersaturation, check_if_output_matches_expected_output) {
    double EPSILON = 0.0001;
    EXPECT_TRUE(fabs(critical_saturation(1.e-6, 273.15) -
                     1.3070837533249676e-05) < EPSILON);
    EXPECT_TRUE(fabs(critical_saturation(1.e-8, 273.15) -
                     0.013070837533249675) < EPSILON);
}

TEST(condensation_solver, check_if_outout_matches_expected_outout) {
    double EPSILON = 0.0001;
    EXPECT_TRUE(fabs(condensation_solver(1.e-6, 1500, 273.15, 0.01, 0, 0.1) -
                     1.1117545636977475e-06) < EPSILON);
    EXPECT_TRUE(fabs(condensation_solver(1.e-6, 1500, 273.15, 0.01, 0, 0.1) -
                     1.1102670529367745e-06) < EPSILON);
    EXPECT_TRUE(fabs(condensation_solver(1.e-8, 1000, 273.15, 0.01, 0, 1) -
                     9.471927239857077e-05) < EPSILON);
}


