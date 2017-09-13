#pragma once

struct Tendencies {
    double dqc;
    double dT;

    Tendencies& operator+=(const Tendencies& rhs){
        this->dqc += rhs.dqc;
        this->dT += rhs.dT;
        return *this;
    }
};
