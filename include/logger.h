#pragma once
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <memory>
#include <fstream>
#include <string>
#include <sstream>
#include <netcdf>
#include <string>
#include "state.h"
#include "superparticle.h"
#include "grid.h"
#include "thermodynamic.h"
#include "analize_sp.h"
#include "analize_state.h"
#include "time_stamp.h"

class Logger {
   public:
    virtual void initialize(const Grid& grid){} 
    virtual void log(const State& state,
                     const std::vector<Superparticle>& superparticles,
                     const int& N_sp) = 0;
    virtual void finalize() const = 0;
    virtual ~Logger(){}
};

class ASCIILogger: public Logger{
    inline void log(const State& state,
                    const std::vector<Superparticle>& superparticles,
                    const int& N_sp) override {
        std::vector<double> qc_sum = calculate_qc_profile(superparticles, state.grid);
        std::vector<double> r_mean = calculate_mean_radius_profile(superparticles, state.grid);
        std::vector<double> r_max = calculate_maximal_radius_profile(superparticles, state.grid);
        std::vector<double> sp_count = count_superparticles(superparticles, state.grid);
        std::vector<int> sp_count_nuc = count_nucleated(superparticles, state.grid);
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
            file << std::setprecision(3) << std::setw(10) << state.grid.getlay(i);
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
                    const int& N_sp)  override {
        std::vector<double> qc_sum = calculate_qc_profile(superparticles, state.grid);
        std::vector<double> r_mean = calculate_mean_radius_profile(superparticles, state.grid);
        std::vector<double> r_max = calculate_maximal_radius_profile(superparticles, state.grid);
        std::vector<double> sp_count = count_superparticles(superparticles, state.grid);
        std::vector<int> sp_count_nuc = count_nucleated(superparticles, state.grid);
        std::vector<double> S = supersaturation_profile(state);

        std::cout << std::endl;
        std::cout << "State at " << state.t << "\n";
        std::cout << "     layer         z         E         p         T       "
                     " qv         S        qc    r_mean     r_max     N_tot     N_nuc\n";
        for (unsigned int i = 0; i < state.layers.size(); ++i) {
            std::cout << std::setprecision(3) << std::setw(10) << i;
            std::cout << std::setprecision(3) << std::setw(10) << state.grid.getlay(i);
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

class NetCDFLogger: public Logger {
    public:
    NetCDFLogger() {
        mkdir(folder.c_str(), S_IRWXU);
        file = time_stamp() + ".nc" ;
        fh = std::make_unique<netCDF::NcFile>(folder+file, netCDF::NcFile::replace);
    }
    virtual void initialize(const Grid& grid){
        n_lay = grid.n_lay;
        netCDF::NcDim layer_dim = fh->addDim("layer", n_lay);
        netCDF::NcVar layer_var = fh->addVar("layer", netCDF::ncDouble, layer_dim);
        auto layers = grid.getlays();
        layer_var.putVar({0}, {n_lay}, layers.data());

        time_dim = fh->addDim("time");
        time_var = fh->addVar("time", netCDF::ncDouble, time_dim);

        qc_var = fh->addVar("qc", netCDF::ncDouble, {time_dim, layer_dim});
        S_var = fh->addVar("S", netCDF::ncDouble, {time_dim, layer_dim});
        r_max_var = fh->addVar("r_max", netCDF::ncDouble, {time_dim, layer_dim});
        r_mean_var = fh->addVar("r_mean", netCDF::ncDouble, {time_dim, layer_dim});
        sp_count_var = fh->addVar("sp_count", netCDF::ncDouble, {time_dim, layer_dim});
        r_std_var = fh->addVar("r_std", netCDF::ncDouble, {time_dim, layer_dim});
        i = 0;
    } 

    inline void log(const State& state,
                    const std::vector<Superparticle>& superparticles,
                    const int& N_sp) override {


        auto qc = calculate_qc_profile(superparticles, state.grid);
        auto S = supersaturation_profile(state);
        auto r_max = calculate_maximal_radius_profile(superparticles, state.grid);
        auto r_mean = calculate_mean_radius_profile(superparticles, state.grid);
        auto sp_count = count_superparticles(superparticles, state.grid);
        auto r_std = calculate_stddev_radius_profile(superparticles, state.grid);

        time_var.putVar({i}, {1}, &state.t);
        qc_var.putVar({i,0}, {1, n_lay}, qc.data());
        S_var.putVar({i,0}, {1, n_lay}, S.data());
        r_max_var.putVar({i,0}, {1, n_lay}, r_max.data());
        r_mean_var.putVar({i,0}, {1, n_lay}, r_mean.data());
        sp_count_var.putVar({i,0}, {1, n_lay}, sp_count.data());
        r_std_var.putVar({i,0}, {1, n_lay}, r_std.data());


        std::cout << "time [s]: " << state.t << std::endl;
        ++i;
    }

    ~NetCDFLogger(){
        //flush();
    }

    inline void finalize() const override {}

    private:
    netCDF::NcDim time_dim;
    netCDF::NcVar time_var;
    netCDF::NcVar qc_var;
    netCDF::NcVar S_var;
    netCDF::NcVar r_max_var;
    netCDF::NcVar r_mean_var;
    netCDF::NcVar sp_count_var;
    netCDF::NcVar r_std_var;
    size_t i;
    size_t n_lay;
    std::string folder = "./data/";
    std::string file;
    std::unique_ptr<netCDF::NcFile> fh;
};

class NetcdfLogger : public Logger {};

inline std::unique_ptr<Logger> createLogger(std::string logger) {
    if(logger == "ascii") {return std::make_unique<ASCIILogger>();}
    if(logger == "std") {return std::make_unique<StdoutLogger>();}
    if(logger == "netcdf") {return std::make_unique<NetCDFLogger>();}
    else {
        return std::make_unique<StdoutLogger>();//sollte eigentlich nicht gehen oder einen std logger aufrufen
    }
}
