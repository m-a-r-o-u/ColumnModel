#pragma once

struct Level{
    double w;
    double p;
    void operator+=(const Level& l){
        w += l.w;
    }
};
