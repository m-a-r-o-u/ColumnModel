#pragma once
#include <vector>
#include <cstddef>
#include <iostream>
#include <iterator>

class LinearField {
   public:
    virtual double operator()(double z, double dz) const = 0;
    virtual LinearField& change(double z, double dz, double dq) = 0;
    virtual LinearField& advect(const LinearField& w, double dt) = 0;
};

class VectorLinearField : public LinearField {
   public:
    inline VectorLinearField(std::vector<double> values, double gridlength)
        : values(values), gridlength(gridlength) {}

    inline VectorLinearField(const VectorLinearField& v_in)
        : values(v_in.values), gridlength(v_in.gridlength) {}

    inline double operator()(double z, double dz) const override {
        return values[z / gridlength];
    }
    inline LinearField& change(double z, double dz, double dv) override {
        values[z / gridlength] += dv * dz / gridlength;
        return *this;
    }
    inline LinearField& advect(const LinearField& w, double dt) override {
        std::vector<double> values_advect(values.begin(), values.end());
        auto scale = dt / gridlength;
        for (size_t i = 1; i < values.size() - 1; ++i) {
            auto a_lo = w(i * gridlength, gridlength);
            auto a_hi = w((i + 1) * gridlength, gridlength);
            if (a_lo < 0) {
                values_advect[i] += scale * a_lo * values[i];
            } else {
                values_advect[i] += scale * a_lo * values[i - 1];
            }
            if (a_hi < 0) {
                values_advect[i] -= scale * a_hi * values[i + 1];
            } else {
                values_advect[i] -= scale * a_hi * values[i];
            }
        }
        std::swap(values, values_advect);
        return *this;
    }
    inline void print_field() {
        std::copy(values.begin(), values.end(),
                  std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;
    }

    inline const double& operator[](size_t i) const { return values[i]; }

    inline auto begin() const {
        return values.begin();
    }
    //auto begin() { return values.begin(); }
    inline std::vector<double>::const_iterator end() const { return values.end(); }

    inline size_t size() const { return values.size(); }

   private:
    std::vector<double> values;
    double gridlength;
};
