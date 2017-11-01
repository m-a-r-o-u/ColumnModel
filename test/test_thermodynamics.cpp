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

TEST(condensation,
     check_if_condensation_works_for_evaporation_S_negative_qc_zero) {
    double EPSILON = 1.e-10;

    double qc = 0.0;
    double N = 1.e8;
    double r_dry = 1.e-6;
    double S = -1;
    double T = 288.;
    double E = 0.;
    double dt = 0.1;
    Tendencies tendencies = condensation(qc, N, r_dry, S, T, E, dt);

    EXPECT_TRUE(std::abs(tendencies.dqc - 0.) < EPSILON);
    EXPECT_TRUE(std::abs(tendencies.dT - 0.) < EPSILON);
}

TEST(
    condensation,
    check_if_condensation_works_for_evaporation_until_r_new_equals_r_dry_S_negative_qc_not_zero) {
    double EPSILON = 1.e-10;

    double qc = 1.e-10;
    double N = 1.e8;
    double r_dry = 1.e-6;
    double S = -1;
    double T = 288.;
    double E = 0.;
    double dt = 0.1;
    Tendencies tendencies = condensation(qc, N, r_dry, S, T, E, dt);

    EXPECT_TRUE(std::abs(tendencies.dqc + qc) < EPSILON);
}

TEST(
    condensation,
    check_if_condensation_works_for_evaporation_until_r_new_equals_r_dry_S_a_little_negative_qc_not_zero) {

    double qc = 1.e-3;
    double N = 1.e8;
    double r_dry = 1.e-6;
    double S = -1.e-10;
    double T = 288.;
    double E = 0.;
    double dt = 0.1;
    Tendencies tendencies = condensation(qc, N, r_dry, S, T, E, dt);

    EXPECT_TRUE(tendencies.dqc > -qc);
    EXPECT_TRUE(tendencies.dqc < 0);
}

TEST(advect_first_order,
     check_advection_of_a_one_with_positive_wind_to_another_cell) {
    std::vector<double> q{1, 0, 0};
    double gridlength = 1;
    double dt = 1;
    std::vector<double> w{1, 1};

    advect_first_order(q.begin(), q.end(), w.begin(), gridlength, dt);
    EXPECT_TRUE(q[0] == 1);
    EXPECT_TRUE(q[1] > 0);
    EXPECT_TRUE(q[1] == 1);
    EXPECT_TRUE(q[2] == 0);
}

TEST(advect_first_order,
     check_advection_of_a_one_with_negative_wind_to_another_cell) {
    std::vector<double> q{1, 0, 1};
    double gridlength = 1;
    double dt = 1;
    std::vector<double> w{-1, -1};

    advect_first_order(q.begin(), q.end(), w.begin(), gridlength, dt);
    EXPECT_TRUE(q[0] == 1);
    EXPECT_TRUE(q[1] > 0);
    EXPECT_TRUE(q[1] == 1);
    EXPECT_TRUE(q[2] == 1);
}

TEST(advect_first_order, check_if_exception_raised_when_cfl_broken) {
    std::vector<double> q{1, 0, 1};
    double gridlength = 1;
    double dt = 2;
    std::vector<double> w{-1, -1};

    EXPECT_ANY_THROW(
        advect_first_order(q.begin(), q.end(), w.begin(), gridlength, dt));
}

TEST(radius_function, check_if_expected_output_is_matched) {
    double EPSILON = 1.e-10;
    double qc = 1.e-3;
    double N = 1.e8;
    double r_min = 1.e-6;
    EXPECT_TRUE(std::abs(radius(qc, N, r_min) - 1.3366912e-5) < EPSILON);
}

TEST(radius_function, check_if_N_is_int) {
    double EPSILON = 1.e-10;
    double qc = 1.e-3;
    int N = 1.e8;
    double r_min = 1.e-6;
    EXPECT_TRUE(std::abs(radius(qc, N, r_min) - 1.3366912e-5) < EPSILON);
}

TEST(radius_function, check_if_qc_is_zero) {
    double EPSILON = 1.e-10;
    double qc = 0;
    int N = 1.e8;
    double r_min = 1.e-6;
    EXPECT_TRUE(std::abs(radius(qc, N, r_min) - r_min) < EPSILON);
}

TEST(first_order_upwind, test_updraft){
    std::vector<double> q{0,1,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;

    first_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 1);
}
TEST(first_order_upwind, test_updraftfromground){
    std::vector<double> q{1,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;

    first_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 1);
    EXPECT_EQ(q[1], 1);
    EXPECT_EQ(q[2], 0);
}
TEST(first_order_upwind, test_updraftfromend){
    std::vector<double> q{0,0,1};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;

    first_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
}
TEST(first_order_upwind, test_dndraft_top){
    std::vector<double> q{0,0,1};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;

    first_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 1);
    EXPECT_EQ(q[2], 1);
}
TEST(first_order_upwind, test_dndraft_middle){
    std::vector<double> q{0,1,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;

    first_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 1);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
}
TEST(first_order_upwind, test_dndraft_ground){
    std::vector<double> q{1,0,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;

    first_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
}
TEST(second_order_upwind, test_updraft_empty){
    std::vector<double> q{0,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
}
TEST(second_order_upwind, test_updraft_ground){
    std::vector<double> q{1,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], 1);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], -1/2.);
}
TEST(second_order_upwind, test_updraft_1lo){
    std::vector<double> q{0,1,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 1);
    EXPECT_EQ(q[2], 2.);
}

TEST(second_order_upwind, test_updraft_cur){
    std::vector<double> q{0,0,1};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], -.5);
}
TEST(second_order_upwind, test_dndraft_2lo){
    std::vector<double> q{0,0,1};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], -.5);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 1);
}
TEST(second_order_upwind, test_dndraft_1lo){
    std::vector<double> q{0,1,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], 2.);
    EXPECT_EQ(q[1], 1);
    EXPECT_EQ(q[2], 0);
}

TEST(second_order_upwind, test_dndraft_cur){
    std::vector<double> q{1,0,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    second_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], -1./2.);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
}

TEST(third_order_upwind, test_updraft_empty){
    std::vector<double> q{0,0,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;

    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
    EXPECT_EQ(q[3], 0);
}
TEST(third_order_upwind, test_updraft_2lo){
    std::vector<double> q{1,0,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 1);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], -1./6.);
    EXPECT_EQ(q[3], 0);
}
TEST(third_order_upwind, test_updraft_1lo){
    std::vector<double> q{0,1,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 1);
    EXPECT_EQ(q[2], 1);
    EXPECT_EQ(q[3], 0);
}
TEST(third_order_upwind, test_updraft_hi){
    std::vector<double> q{0,0,0,1};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], -1/3.);
    EXPECT_EQ(q[3], 1);
}
TEST(third_order_upwind, test_dndraft_empty){
    std::vector<double> q{0,0,0,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
    EXPECT_EQ(q[3], 0);
}
TEST(third_order_upwind, test_dndraft_2lo){
    std::vector<double> q{0,0,0,1};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], -1/6.);
    EXPECT_EQ(q[2], 0);
    EXPECT_EQ(q[3], 1);
}
TEST(third_order_upwind, test_dndraft_1lo){
    std::vector<double> q{0,0,1,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 1);
    EXPECT_EQ(q[2], 1);
    EXPECT_EQ(q[3], 0);
}
TEST(third_order_upwind, test_dndraft_hi){
    std::vector<double> q{1,0,0,0};
    std::vector<double> w{-1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 1);
    EXPECT_EQ(q[1], -1/3.);
    EXPECT_EQ(q[2], 0);
    EXPECT_EQ(q[3], 0);
}

TEST(sixth_order_wickerskamarock, test_updraft_empty){
    std::vector<double> q{0,0,0,0,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    third_order_upwind(q.begin(), q.end(), w.begin(), dl, dt);
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
    EXPECT_EQ(q[3], 0);
}

TEST(sixth_order_wickerskamarock, test_updraft_ground){
    std::vector<double> q{1,0,0,0,0,0};
    std::vector<double> w{1};
    double dl = 1;
    double dt = 1;
    sixth_order_wickerskamarock(q.begin(), q.end(), w.begin(), dl, dt);
    for ( auto qi: q){
        std::cout << qi << std::endl;
    }
    EXPECT_EQ(q[0], 1);
    EXPECT_EQ(q[1], 0);
    EXPECT_EQ(q[2], 0);
    EXPECT_EQ(q[3], -1/60.);
    EXPECT_EQ(q[4], 0);
    EXPECT_EQ(q[5], 0);
}
