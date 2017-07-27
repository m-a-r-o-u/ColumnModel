#pragma once

struct Layer{
    double T;
    double p;
    double qv;
    double E;
    void operator+=(const Layer& l){
        T += l.T;
        p += l.p;
        qv += l.qv;
        E += l.E;
    }
};
