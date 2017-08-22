#pragma once
#include <iomanip>
#include <memory>
#include <fstream>
#include <string>
#include <sstream>
#include "state.h"
#include "superparticle.h"
#include "grid.h"
#include "thermodynamic.h"
#include "analize_sp.h"
#include "analize_state.h"

class Logger {
   public:
    virtual void log(const State& state,
                     const std::vector<Superparticle>& superparticles,
                     const Grid& grid, const int& N_sp) const = 0;
    virtual void finalize() const = 0;
};

class ASCIILogger: public Logger{
    inline void log(const State& state,
                    const std::vector<Superparticle>& superparticles,
                    const Grid& grid, const int& N_sp) const override {
        std::vector<double> qc_sum = calculate_qc_profile(superparticles, grid);
        std::vector<double> r_mean = calculate_mean_radius_profile(superparticles, grid);
        std::vector<double> r_max = calculate_maximal_radius_profile(superparticles, grid);
        std::vector<double> sp_count = count_superparticles(superparticles, grid);
        std::vector<double> sp_count_nuc = count_nucleated(superparticles, grid);
        std::vector<double> S = supersaturation_profile(state);


        std::fstream file;
        std::stringstream file_name_stream;

        file_name_stream << "./out";
//        file_name_stream << int(superparticles[0].N) << ".";
        file_name_stream << '.' << std::setfill('0') << std::setw(5) << int(state.t);
        file_name_stream << '.' << std::setfill('0') << std::setw(5) << int(N_sp);

        std::string file_name;
        file_name = file_name_stream.str();

        file.open(file_name, std::ios::out);

        for (unsigned int i = 0; i < state.layers.size(); ++i) {
            file << std::setprecision(3) << std::setw(10) << i;
            file << std::setprecision(3) << std::setw(10) << grid.getlay(i);
            file << std::setprecision(3) << std::setw(10) << state.layers[i].E;
            file << std::setprecision(3) << std::setw(10) << state.layers[i].p;
            file << std::setprecision(3) << std::setw(10) << state.layers[i].T;
            file << std::setprecision(3) << std::setw(10) << state.layers[i].qv;
            file << std::setprecision(3) << std::setw(10) << S[i];
            file << std::setprecision(3) << std::setw(10) << qc_sum[i];
            file << std::setprecision(3) << std::setw(10) << r_mean[i];
            file << std::setprecision(3) << std::setw(10) << r_max[i];
            file << std::setprecision(3) << std::setw(10) << sp_count_nuc[i];
            file << "\n";
        }
        file.close();
    }
    inline void finalize() const override {}

};

class StdoutLogger : public Logger {
   public:
    inline void log(const State& state,
                    const std::vector<Superparticle>& superparticles,
                    const Grid& grid, const int& N_sp) const override {
        std::vector<double> qc_sum = calculate_qc_profile(superparticles, grid);
        std::vector<double> r_mean = calculate_mean_radius_profile(superparticles, grid);
        std::vector<double> r_max = calculate_maximal_radius_profile(superparticles, grid);
        std::vector<double> sp_count = count_superparticles(superparticles, grid);
        std::vector<double> sp_count_nuc = count_nucleated(superparticles, grid);
        std::vector<double> S = supersaturation_profile(state);

        std::cout << std::endl;
        std::cout << "State at " << state.t << "\n";
        std::cout << "     layer         z         E         p         T       "
                     " qv         S        qc    r_mean     r_max     N_tot     N_nuc\n";
        for (unsigned int i = 0; i < state.layers.size(); ++i) {
            std::cout << std::setprecision(3) << std::setw(10) << i;
            std::cout << std::setprecision(3) << std::setw(10) << grid.getlay(i);
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].E;
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].p;
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].T;
            std::cout << std::setprecision(3) << std::setw(10) << state.layers[i].qv;
            std::cout << std::setprecision(3) << std::setw(10) << S[i];
            std::cout << std::setprecision(3) << std::setw(10) << qc_sum[i];
            std::cout << std::setprecision(3) << std::setw(10) << r_mean[i];
            std::cout << std::setprecision(3) << std::setw(10) << r_max[i];
            std::cout << std::setprecision(3) << std::setw(10) << sp_count[i];
            std::cout << std::setprecision(3) << std::setw(10) << sp_count_nuc[i];
            std::cout << "\n";
        }
        std::cout << std::endl;
    }
    inline void finalize() const override {}
};

class NetcdfLogger : public Logger {};

inline std::unique_ptr<Logger> createLogger(std::string logger) {
    if(logger == "ascii") {return std::make_unique<ASCIILogger>();}
    else {return std::make_unique<StdoutLogger>();}

}
