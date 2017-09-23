#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "file_utils.h"
#include "twomey_utils.h"

void load_data(std::vector<double>& n, std::vector<double>& s);

void set_n(std::vector<double>& n_out, const std::vector<double>& n, double N);

double linear_interpolate(double x1, double y1, double x2, double y2, double x);

std::vector<double> calculate_stable(const std::vector<double>& n, const std::vector<double>& s, const std::vector<double>& nx);

std::vector<double> nstable(int N, int& Nmulti);

template <typename T, typename U, typename V>
std::vector<double> arange(T start, U stop, V step){
    std::vector<double> out;
    int i = 0;
    double item = start + i * step;
    if (step > 0){
        while(item < stop){
            out.push_back(item);
            ++i;
            item = start + i * step;
        }
    }
    return out;
}

