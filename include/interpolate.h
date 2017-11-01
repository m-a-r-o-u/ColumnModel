#pragma once
#include <iostream>

inline double linear_interpolate(double x1, double y1, double x2, double y2, double x){
    if (x1 > x2) {
        std::cout << "linear interpolation: points may be in the wrong order: x1 "<<  x1 << " x2 " << x2 << std::endl;
    }
    double a = ( y2 - y1 ) / ( x2 - x1 );
    double b = y2 - a * x2;
    return a * x + b;
}

inline double bi_linear_interpolate(double x11, double x12, double v11, double v12, double x21, double x22, double v21, double v22, double x, double y){
    if (x11 > x21) {
        std::cout << "linear interpolation: points may be in the wrong order: x11 "<<  x11 << " x21 " << x21 << std::endl;
    }
    if (x12 > x22) {
        std::cout << "linear interpolation: points may be in the wrong order: x12 "<<  x12 << " x22 " << x22 << std::endl;
    }
    double vinter1 =  linear_interpolate(x11, v11, x21, v21, x);
    double vinter2 =  linear_interpolate(x11, v12, x21, v22, x);
    return linear_interpolate(x12, vinter1, x22, vinter2, y);
}
