#pragma once
#include "linearfield.h"

struct State {
    double t;
    VectorLinearField w;
    VectorLinearField p;
    VectorLinearField T;
    VectorLinearField E;
    VectorLinearField qv;
};
