#pragma once
#include <iterator>

template <typename BaseIt, typename Element>
class MemberIterator {
   public:
    typedef
        typename std::iterator_traits<BaseIt>::difference_type difference_type;
    typedef Element value_type;
    typedef Element* pointer;
    typedef Element& reference;
    typedef std::forward_iterator_tag iterator_category;
    typedef typename std::iterator_traits<BaseIt>::value_type StructT;

    MemberIterator(BaseIt base_it, Element StructT::*member_ptr)
        : base_it(base_it), member_ptr(member_ptr){};

    MemberIterator<BaseIt, Element>& operator++() {
        ++base_it;
        return *this;
    }
    Element& operator*() { return (*base_it).*member_ptr; }
    bool operator!=(const MemberIterator<BaseIt, Element>& other) {
        return base_it != other.base_it;
    }
    int operator-(const MemberIterator<BaseIt, Element>& other) {
        return base_it - other.base_it;
    }

   private:
    BaseIt base_it;
    Element StructT::*member_ptr;
};

template <typename BaseIt, typename Element>
MemberIterator<BaseIt, Element> member_iterator(
    BaseIt base_it,
    Element std::iterator_traits<BaseIt>::value_type::*member_ptr) {
    return {base_it, member_ptr};
}
