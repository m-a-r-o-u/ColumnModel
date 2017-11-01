#pragma once
#include <ostream>

struct Tendencies {
    double dqc;
    double dT;

    Tendencies& operator+=(const Tendencies& rhs){
        this->dqc += rhs.dqc;
        this->dT += rhs.dT;
        return *this;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Tendencies& rhs){
    os << rhs.dqc << " " << rhs.dT;
    return os;
}
