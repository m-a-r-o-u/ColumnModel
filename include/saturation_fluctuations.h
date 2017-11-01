#pragma once
#include <limits>
#include <random>
#include "superparticle.h"
#include "constants.h"
#include "tau_relax.h"

inline double turbulent_kinetic_energy(const double& l, const double& epsilon) {
    double c = 0.845;
    return std::pow(l * epsilon / c, 2. / 3.);
}

inline double w_standart(const double& e) { return std::sqrt(2. / 3. * e); }

inline double integral_timescale(const double& l, const double& tke) {
    const double c = 1.5;
    return l / std::pow(2. * PI, 1./3.) * std::sqrt(c / tke);
}

template <typename G>
double ornstein_uhlenbeck_process(G& gen,const double& w, const double& dt,
                                         const double& tau,
                                         const double& w_std) {
    std::normal_distribution<> d(0., 1.);
    return w * std::exp(-dt / tau) +
           std::sqrt(1 - std::exp(-2 * dt / tau)) * w_std * d(gen);
}

inline double saturation_fluctuations(const double& w_prime, const double& dt,
                                      const double& tau_r,
                                      const double& S_prime) {
    double a1 = 3.e-4;
    return S_prime + dt * (a1 * w_prime - S_prime / tau_r);
}

class FluctuationSolver {
   public:
    virtual void refresh(const std::vector<Superparticle>& sp) = 0;
    virtual double getFluctuation(Superparticle& s, const double& dt) = 0;
};

template <typename G>
class MarkovFluctuationSolver : public FluctuationSolver {
   public:
    MarkovFluctuationSolver(G& gen, const double& epsilon, double l, const Grid& grid)
        : epsilon(epsilon), l(l), gen(gen), tau_relax(grid) {}
    void refresh(const std::vector<Superparticle>& sp) override;
    double getFluctuation(Superparticle& s, const double& dt) override;

   private:
    const double epsilon;
    const double l;
    G& gen;
    TauRelax tau_relax;
};

template <typename G>
void MarkovFluctuationSolver<G>::refresh(
    const std::vector<Superparticle>& sp) {
    tau_relax.refresh(sp);
}

template <typename G>
double MarkovFluctuationSolver<G>::getFluctuation(Superparticle& s,
                                                      const double& dt) {
    double tke = turbulent_kinetic_energy(l, epsilon);
    double tau = integral_timescale(l, tke);
    double w_std = w_standart(tke);
    s.w_prime = ornstein_uhlenbeck_process(gen, s.w_prime, dt, tau, w_std);
    double tau_r = tau_relax(s.z);
    s.S_prime = saturation_fluctuations(s.w_prime, dt, tau_r, s.S_prime);
    return s.S_prime;
}

class NoFluctuationSolver : public FluctuationSolver {
   public:
    NoFluctuationSolver(){}
    void refresh(const std::vector<Superparticle>& sp) override {}
    double getFluctuation(Superparticle& s, const double& dt) override {return 0.;}
};

template <typename G>
std::unique_ptr<FluctuationSolver> mkFS(G& gen, 
        const std::string& type, const double& epsilon, double l, const Grid& grid) {
    if (type == "markov") {
        return std::make_unique<MarkovFluctuationSolver<G>>(gen, epsilon, l, grid);
    }
    else {
        return std::make_unique<NoFluctuationSolver>();
    }
//    else{
//        //throw no fluctuation exceptioon
//    }
}
