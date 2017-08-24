#pragma once
#include <algorithm>
#include <vector>
#include "thermodynamic.h"
#include "state.h"

inline std::vector<double> supersaturation_profile(const State& state) {
    std::vector<double> res(state.layers.size(), 0.);
    std::transform(state.layers.begin(), state.layers.end(), res.begin(),
                   [](Layer l) { return saturation(l.T, l.p, l.qv); });
    return res;
}
