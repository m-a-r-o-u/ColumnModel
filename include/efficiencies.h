#pragma once
#include "interpolate.h"
#include "twomey_utils.h"

class Efficiencies {
   public:
    Efficiencies() {}
    double collision_efficiency(double R, double rR) const {
        unsigned int iR = R * 0.1;
        if (iR >= sizeof(Rref_remap) / sizeof(unsigned int)) {
            iR = sizeof(Rref_remap) / sizeof(unsigned int) - 1;
        }
        iR = Rref_remap[iR];
        // auto iR = left_index_min_zero_max_smallerlast(Rref, R);
        unsigned int irR = (rR * 20) - 1;
        if (irR >= 19) {
            irR = 18;
        }
        // auto irR = left_index_min_zero_max_smallerlast(rRref, rR);
        return bi_linear_interpolate(
            rRref[irR], Rref[iR], efficiencies[iR][irR],
            efficiencies[iR + 1][irR], rRref[irR + 1], Rref[iR + 1],
            efficiencies[iR][irR + 1], efficiencies[iR + 1][irR + 1], rR, R);
    }

   private:
    unsigned char Rref_remap[31] = {
        0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8,
        8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9};
    // std::vector<double> Rref;
    // std::vector<double> rRref;
    // std::vector<std::vector<double>> efficiencies;

    double Rref[11] = {10, 20,  30,  40,  50, 60,
                                        70, 100, 150, 200, 300};
    double rRref[20] = {
        0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
        0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
    double efficiencies[11][20] = {
        {0.0001, 0.0001, 0.0001, 0.014, 0.017, 0.019, 0.022,
         0.027,  0.030,  0.033,  0.035, 0.037, 0.038, 0.038,
         0.037,  0.036,  0.035,  0.032, 0.029, 0.027},
        {0.0001, 0.0001, 0.005, 0.016, 0.022, 0.03,  0.043,
         0.052,  0.064,  0.072, 0.079, 0.082, 0.080, 0.076,
         0.067,  0.057,  0.048, 0.040, 0.033, 0.027},
        {0.0001, 0.002, 0.02, 0.04, 0.085, 0.17, 0.27, 0.40, 0.50, 0.55,
         0.58,   0.59,  0.58, 0.54, 0.51,  0.49, 0.47, 0.45, 0.47, 0.52},
        {0.001, 0.07, 0.28, 0.50, 0.62, 0.68, 0.74, 0.78, 0.80, 0.80,
         0.80,  0.78, 0.77, 0.76, 0.77, 0.77, 0.78, 0.79, 0.95, 1.40},
        {0.005, 0.40, 0.60, 0.70, 0.78, 0.83, 0.86, 0.88, 0.90, 0.90,
         0.90,  0.90, 0.89, 0.88, 0.88, 0.89, 0.92, 1.01, 1.30, 2.30},
        {0.05, 0.43, 0.64, 0.77, 0.84, 0.87, 0.89, 0.90, 0.91, 0.91,
         0.91, 0.91, 0.91, 0.92, 0.93, 0.95, 1.00, 1.03, 1.70, 3.00},
        {0.20, 0.58, 0.75, 0.84, 0.88, 0.90, 0.92, 0.94, 0.95, 0.95,
         0.95, 0.95, 0.95, 0.95, 0.97, 1.00, 1.02, 1.04, 2.30, 4.00},
        {0.50, 0.79, 0.91, 0.95, 0.95, 1.00, 1.00, 1.00, 1.00, 1.00,
         1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00},
        {0.77, 0.93, 0.97, 0.97, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
         1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00},
        {0.87, 0.96, 0.98, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
         1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00},
        {0.97, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
         1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00}};
};