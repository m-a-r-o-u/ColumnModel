#pragma once
#include <algorithm>
#include <vector>

template <typename IT, typename T>
int lower_bound_index(IT first, IT last, T value) {
    return std::distance(first, std::lower_bound(first, last, value));
}

template <typename IT, typename T>
int upper_bound_index(IT first, IT last, T value) {
    return std::distance(first, std::upper_bound(first, last, value));
}

template <typename C>
std::vector<int> indexes(const C& Slvl, const C& Sprf) {
    std::vector<int> Nprf;
    Nprf.reserve(Sprf.size());
    for (auto s : Sprf) {
        Nprf.push_back(upper_bound_index(Slvl.begin(), Slvl.end(), s));
    }
    return Nprf;
}

template <typename C, typename T>
int left_index(const C& n, T nx){
    int i = lower_bound_index(n.begin(), n.end(), nx) - 1;
    return i;
}

template <typename C, typename T>
int right_index(const C& n, T nx){
    return lower_bound_index(n.begin(), n.end(), nx);
}

template <typename C, typename T>
int left_index_min_zero( const C& n, const T nx){
    int li = left_index(n, nx);
    if (li == -1){
        li = 0;
    }
    return li;
}

template <typename C, typename T>
int left_index_min_zero_max_smallerlast( const C& n, const T nx){
    int li = left_index(n, nx);
    if (li == -1){
        li = 0;
    }
    if (li == int(n.size()) - 1 ){
        --li;
    }
    return li;
}

