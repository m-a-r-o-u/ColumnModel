#pragma once
#include <iomanip>
#include <iostream>

struct BackgroundLevelAfglus {
    double z;
    double p;
    double T;
    double air;
    double o3;
    double o2;
    double h2o;
    double co2;
    double no2;
};

std::ostream& operator<<(std::ostream& os, const BackgroundLevelAfglus& bl) {
    int w = 10;
    int p = 3;
    return os << std::setprecision(p) << std::setw(w) << bl.z
              << std::setprecision(p) << std::setw(w) << bl.p
              << std::setprecision(p) << std::setw(w) << bl.T
              << std::setprecision(p) << std::setw(w) << bl.air
              << std::setprecision(p) << std::setw(w) << bl.o3
              << std::setprecision(p) << std::setw(w) << bl.o2
              << std::setprecision(p) << std::setw(w) << bl.h2o
              << std::setprecision(p) << std::setw(w) << bl.co2
              << std::setprecision(p) << std::setw(w) << bl.no2;
}
std::istream& operator>>(std::istream& is, BackgroundLevelAfglus& bl) {
    double dummy;  // get rid by using ignore or a similar std function
    return is >> bl.z >> bl.p >> bl.T >> bl.air >> bl.o3 >> bl.o2 >> bl.h2o >>
           bl.co2 >> bl.no2;
}
