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
    typedef typename std::iterator_traits<BaseIt>::iterator_category iterator_category;
    typedef typename std::iterator_traits<BaseIt>::value_type StructT;

    MemberIterator(BaseIt base_it, Element StructT::*member_ptr)
        : base_it(base_it), member_ptr(member_ptr){};

    MemberIterator<BaseIt, Element>& operator++() {
        ++base_it;
        return *this;
    }
    MemberIterator<BaseIt, Element> operator++(int) {
        auto that = *this;
        ++base_it;
        return that;
    }
    MemberIterator<BaseIt, Element>& operator--() {
        --base_it;
        return *this;
    }
    MemberIterator<BaseIt, Element> operator--(int) {
        auto that = *this;
        --base_it;
        return that;
    }
    Element& operator*() { return (*base_it).*member_ptr; }
    Element* operator->() {
        return &((*base_it).*member_ptr);
    }
    Element& operator[](difference_type n){
        return base_it[n].*member_ptr;
    }
    bool operator<(const MemberIterator<BaseIt, Element>& other) const {
        return base_it < other.base_it;
    }
    bool operator<=(const MemberIterator<BaseIt, Element>& other) const {
        return base_it <= other.base_it;
    }
    bool operator>(const MemberIterator<BaseIt, Element>& other) const {
        return base_it > other.base_it;
    }
    bool operator>=(const MemberIterator<BaseIt, Element>& other) const {
        return base_it >= other.base_it;
    }
    bool operator!=(const MemberIterator<BaseIt, Element>& other) const {
        return base_it != other.base_it;
    }
    bool operator==(const MemberIterator<BaseIt, Element>& other) const {
        return base_it == other.base_it;
    }
    MemberIterator<BaseIt, Element> operator+(difference_type n) const {
        return {base_it + n, member_ptr};
    }
    MemberIterator<BaseIt, Element>& operator+=(difference_type n) {
        base_it += n;
        return *this;
    }
    MemberIterator<BaseIt, Element>& operator-=(difference_type n) {
        base_it -= n;
        return *this;
    }
    friend MemberIterator<BaseIt, Element> operator+(difference_type n, const MemberIterator<BaseIt, Element>& rhs) {
        return rhs + n;
    }
    MemberIterator<BaseIt, Element> operator-(difference_type n) const {
        return {base_it - n, member_ptr};
    }
    difference_type operator-(const MemberIterator<BaseIt, Element>& other) {
        return base_it - other.base_it;
    }
   private:
    BaseIt base_it;
    Element StructT::* member_ptr;
};

template <typename BaseIt, typename Element>
MemberIterator<BaseIt, Element> member_iterator(
    BaseIt base_it,
    Element std::iterator_traits<BaseIt>::value_type::*member_ptr) {
    return {base_it, member_ptr};
}
