#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

struct Grid {
    Grid(double h, double l) : height(h), length(l){
        for (int i = 0; i < n_lay; ++i) {
            lays.push_back((i + 0.5) * length);
        }
        for (int i = 0; i < n_lvl; ++i) {
            lvls.push_back(i * length);
        }
    }
    const double height;
    const double length;
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
    const std::vector<double>& getlays() const {
        return lays;
    }
    const std::vector<double>& getlvls() const {
        return lvls;
    }
    private:
    std::vector<double> lays;
    std::vector<double> lvls;
};
