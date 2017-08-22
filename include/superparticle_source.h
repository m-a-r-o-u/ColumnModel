#pragma once
#include "superparticle.h"

template <typename OutputIterator>
class SuperParticleSource {
   public:
    virtual ~SuperParticleSource() {}
    virtual void generateParticles(OutputIterator it, State& state, double dt,
                                   const std::vector<Superparticle>& sp) = 0;
};

template <typename D, typename G, typename OutputIterator>
class SuperParticleSourceConstHeight
    : public SuperParticleSource<OutputIterator> {
   public:
    SuperParticleSourceConstHeight(double z_insert, double rate, int N, D d,
                                   G g)
        : z_insert(z_insert), rate(rate), N(N), d(d), g(g){};
    void generateParticles(OutputIterator it, State& state,
                           double dt, const std::vector<Superparticle>& sp) override {
        for (int i = 0; i < dt * rate; ++i) {
            *it = Superparticle{0, z_insert, d(g), N, false};
            ++it;
        }
    }

   private:
    double z_insert;
    double rate;
    int N;
    D d;
    G g;
};

template <typename OIt, typename D, typename G>
std::unique_ptr<SuperParticleSourceConstHeight<D, G, OIt>> mkSPSCH(
    double z_insert, double rate, int N, D d, G g) {
    return std::make_unique<SuperParticleSourceConstHeight<D, G, OIt>>(
        z_insert, rate, N, d, g);
}
