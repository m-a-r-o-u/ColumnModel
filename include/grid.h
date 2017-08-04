#pragma once
#include <stdexcept>
#include <cmath>
#include <vector>

struct Grid {
    const double height;
    const double length;
    const double z0;
    const int n_lay = std::ceil(height / length);
    const int n_lvl = n_lay + 1;

    double getlvl(const int i) const {
        if (i < n_lvl && i >= 0) {
            return length * i;
        } else {
            throw std::out_of_range(
                "the lvl index is larger then the number of levels");
        }
    }

    double getlay(const int i) const {
        if (i < n_lay && i >= 0) {
            return length * (i + 0.5);
        } else {
            throw std::out_of_range(
                "the lay index is larger then the number of layers");
        }
    }
    std::vector<double> getlays() const {
        std::vector<double> z;
        for (int i = 0; i < n_lay; ++i) {
            z.push_back((i + 0.5) * length);
        }
        return z;
    }
    std::vector<double> getlvls() const {
        std::vector<double> z;
        for (int i = 0; i < n_lvl; ++i) {
            z.push_back(i * length);
        }
        return z;
    }
};
