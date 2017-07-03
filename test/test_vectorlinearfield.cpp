#include "gtest/gtest.h"
#include "linearfield.h"

TEST(vector_linear_field, copy_is_deep) {
    VectorLinearField vec1({0}, 1);
    VectorLinearField vec2(vec1);
    vec2.change(0.5, 1, 1);

    EXPECT_EQ(vec1(0.5, 1), 0);
    EXPECT_EQ(vec2(0.5, 1), 1);
}
