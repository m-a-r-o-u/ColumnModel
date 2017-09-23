#pragma once
#include <vector>
#include "grid.h"
#include "superparticle.h"

class TauRelax {
   public:
    TauRelax(const Grid& grid)
        : grid(grid) {
    }

    void refresh(const std::vector<Superparticle>& sp);
    inline double operator()(double z) const;

   private:
    std::vector<double> tau_relax;
    const Grid& grid;
};

inline double TauRelax::operator()(double z) const { return tau_relax[grid.getlayindex(z)]; }

