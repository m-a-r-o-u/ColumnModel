#include "tau_relax.h"
#include <algorithm>
#include <limits>

void TauRelax::refresh(const std::vector<Superparticle>& sp) {
    tau_relax.resize(grid.n_lay);
    double a2 = 2.8e-4;
    std::vector<double> one_over_tau(grid.n_lay, 0.);
    for (auto s : sp) {
        int index = grid.getlayindex(s.z);
        one_over_tau[index] += s.radius() * s.N;
    }
    std::transform(one_over_tau.begin(), one_over_tau.end(), tau_relax.begin(),
    [a2](double x) {
        if (x > 0.) {
            return 1. / (x * a2);
        } else {
            return std::numeric_limits<double>::infinity();
        }
                   });
}
