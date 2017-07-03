#pragma once
#include "cmath"

constexpr double PI = std::atan(1)*4; ///< PI value
constexpr double RHO_AIR = 1; ///< [kg m-3] density of air mixture APPROXIMATION
constexpr double RHO_H2O = 1000; ///< [kg m-3] density of water at T=288 [k] and p=100000 [Pa]
constexpr double ES0 = 611.2; ///< [Pa] saturation pressure over flat surface at T=273.15[k]
constexpr double T0 = 273.15; ///< [K] freezing temperature of water
constexpr double H_LAT = 2.257e6; ///< [J kg-1] heat of vaporization
constexpr double R_V = 461.401; ///< [J kg-1 K-1] specific gas constant of water vapor
constexpr double R_G = 287.102; ///< [J kg-1 K-1] ideal gas constant of air
constexpr double K = 2.43e-2; ///< [W m-1 K-1] thermal conductivity of air
constexpr double D = 2.82e-5; ///< [m2 s-1] diffusion constant of water vapor in air
constexpr double C_P = 1003.5; ///< [J kg-1 K-1] specific heat at constant pressure
constexpr double GAMMA = 72.7e-3; ///< [N m-1] surface tension of water at 293 [K]
constexpr double M_MOL_H2O = 18.e-3; ///< [kg mol-1] molecular mass of water
constexpr double M_MOL_S = 58.4e-3; ///< [kg mol-1] molecular mass of NaCl
constexpr double RHO_S = 2.16e3; ///< [kg m-3] density of NaCL at 298 [K]
constexpr double ETA_AIR = 17.1e-6; ///< [Pa s'] dynamic viscosity of air at 273 [K]
constexpr double G = 9.81; ///< [kg m2 s-2] gravity constant
constexpr double LAPSE_RATE_A = G / C_P; ///< [K m-1] adiabatic lapse rate
