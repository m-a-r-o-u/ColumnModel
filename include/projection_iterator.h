#pragma once
#include <iterator>
#include <utility>
#include <type_traits>

template <typename BaseIt, typename F>
class ProjectionIterator {
   public:
    typedef
        typename std::iterator_traits<BaseIt>::difference_type difference_type;
    typedef typename std::iterator_traits<BaseIt>::value_type StructT;
    typedef decltype (std::declval<F>()(std::declval<StructT>())) return_type;
    typedef typename std::remove_reference<return_type>::type value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef typename std::iterator_traits<BaseIt>::iterator_category iterator_category;

    ProjectionIterator(BaseIt base_it, F f)
        : base_it(base_it), f(f){};

    ProjectionIterator<BaseIt, F>& operator++() {
        ++base_it;
        return *this;
    }
    ProjectionIterator<BaseIt, F> operator++(int) {
        auto that = *this;
        ++base_it;
        return that;
    }
    ProjectionIterator<BaseIt, F>& operator--() {
        --base_it;
        return *this;
    }
    ProjectionIterator<BaseIt, F> operator--(int) {
        auto that = *this;
        --base_it;
        return that;
    }
    return_type operator*() { return f(*base_it); }
#pointer operator->() {
        return &f(*base_it);
    }
    return_type operator[](difference_type n){
        return f(base_it[n]);
    }
    bool operator<(const ProjectionIterator<BaseIt, F>& other) const {
        return base_it < other.base_it;
    }
    bool operator<=(const ProjectionIterator<BaseIt, F>& other) const {
        return base_it <= other.base_it;
    }
    bool operator>(const ProjectionIterator<BaseIt, F>& other) const {
        return base_it > other.base_it;
    }
    bool operator>=(const ProjectionIterator<BaseIt, F>& other) const {
        return base_it >= other.base_it;
    }
    bool operator!=(const ProjectionIterator<BaseIt, F>& other) const {
        return base_it != other.base_it;
    }
    bool operator==(const ProjectionIterator<BaseIt, F>& other) const {
        return base_it == other.base_it;
    }
    ProjectionIterator<BaseIt, F> operator+(difference_type n) const {
        return {base_it + n, f};
    }
    ProjectionIterator<BaseIt, F>& operator+=(difference_type n) {
        base_it += n;
        return *this;
    }
    ProjectionIterator<BaseIt, F>& operator-=(difference_type n) {
        base_it -= n;
        return *this;
    }
    friend ProjectionIterator<BaseIt, F> operator+(difference_type n, const ProjectionIterator<BaseIt, F>& rhs) {
        return rhs + n;
    }
    ProjectionIterator<BaseIt, F> operator-(difference_type n) const {
        return {base_it - n, f};
    }
    difference_type operator-(const ProjectionIterator<BaseIt, F>& other) {
        return base_it - other.base_it;
    }
   private:
    BaseIt base_it;
    F f;
};

template <typename BaseIt, typename F>
ProjectionIterator<BaseIt, F> projection_iterator(
    BaseIt base_it,
    F f) {
    return {base_it, f};
}
