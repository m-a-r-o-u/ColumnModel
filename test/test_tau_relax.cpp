#include "gtest/gtest.h"
#include "tau_relax.h"
#include "grid.h"
#include "superparticle.h"

TEST(tau_relax, test_refresh){
    Grid  grid{300., 100.};
    TauRelax relax(grid);

    std::vector<Superparticle> sp;
    sp.push_back({0.00001, 50, 1.e-6, 100000000, true});

    relax.refresh(sp);
    EXPECT_TRUE(true);
}
