#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "member_iterator.h"
#include "constants.h"

template <typename ClassT, typename IFS, typename OI>
void readin_atm(IFS& ifs, OI oi) {
    if (!ifs.good()) {
        std::cout << "ifstream not good, check file name" << std::endl;
    }

    while (ifs.good()) {
        std::string line;
        std::getline(ifs, line);
        if (line.empty() || line.front() == '#') {
            continue;
        }
        std::stringstream linestream(line);
        ClassT ct;
        linestream >> ct;

        *oi = ct;
        ++oi;
    }
}

template <typename Container, typename Element>
std::vector<Element> get_quantity(Container& c,
                                  Element Container::value_type::*mem_ptr) {
    return std::vector<Element>(member_iterator(c.begin(), mem_ptr),
                                member_iterator(c.end(), mem_ptr));
}

template <typename Container>
std::vector<double> pairwise_mean(Container& c) {
    std::vector<double> out;
    out.reserve(c.size() - 1);
    std::transform(c.begin(), c.end() - 1, c.begin() + 1,
                   std::back_inserter(out),
                   [](auto a, auto b) { return (a + b) / 2.; });
    return out;
}

template <typename Container>
std::vector<double> divide_elementwise(Container& nom, Container& denom) {
    std::vector<double> out;
    out.reserve(nom.size());
    std::transform(nom.begin(), nom.end(), denom.begin(),
                   std::back_inserter(out),
                   [](double a, double b) { return a / b; });
    return out;
}

template <typename Container, typename NumberT>
std::vector<double> multiply_all(Container& c, NumberT factor) {
    std::vector<double> out;
    out.reserve(c.size());
    std::transform(c.begin(), c.end(), std::back_inserter(out),
                   [factor](double a) { return a * factor; });
    return out;
}

template <typename Container, typename NumberT>
std::vector<double> from_number_to_volume_ratio(Container& c, NumberT m_mol, NumberT rho) {
    double factor = m_mol / M_MOL_AIR * RHO_AIR / rho;
    std::vector<double> out;
    out.reserve(c.size());
    std::transform(c.begin(), c.end(), std::back_inserter(out),
                   [factor](double a) { return a * factor; });
    return out;
}
