#include "gtest/gtest.h"
#include "member_iterator.h"
#include <vector>
#include <list>

struct Dummy {
    int x;
    int y;
};

TEST(member_iterator, test_lt_op){
    std::vector<Dummy> d{{1,2},{3,4}};
    EXPECT_LT( member_iterator(d.begin(), &Dummy::x), member_iterator(d.end(), &Dummy::x));
}
TEST(member_iterator, test_gt_op){
    std::vector<Dummy> d{{1,2},{3,4}};
    EXPECT_GT( member_iterator(d.end(), &Dummy::x), member_iterator(d.begin(), &Dummy::x));
}
TEST(member_iterator, test_eq_op){
    std::vector<Dummy> d{{1,2},{3,4}};
    EXPECT_EQ( member_iterator(d.begin(), &Dummy::x), member_iterator(d.begin(), &Dummy::x));
}
TEST(member_iterator, test_element_access_op){
    std::vector<Dummy> d{{1,2},{3,4}};
    EXPECT_EQ(member_iterator(d.begin(), &Dummy::x)[1], 3);
}
TEST(member_iterator, test_add_difftype_op){
    std::vector<Dummy> d{{1,2},{3,4}};
    auto tmp = member_iterator(d.begin(), &Dummy::x);
    EXPECT_EQ(tmp+=1, member_iterator(d.begin()+1, &Dummy::x));
    EXPECT_EQ(tmp, member_iterator(d.begin()+1, &Dummy::x));
}
TEST(member_iterator, test_list_deref){
    std::list<Dummy> d{{1,2},{3,4}};
    EXPECT_EQ(*member_iterator(d.begin(), &Dummy::x), 1);
}
