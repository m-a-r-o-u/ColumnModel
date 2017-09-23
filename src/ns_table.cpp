#include "ns_table.h"

void load_data(std::vector<double>& n, std::vector<double>& s){
    auto infile = loadfile("/home/m/Mares.Barekzai/phd/projects/column_model/data/ns_data.txt");
    std::string line;
    while(std::getline(infile, line)){
        if(line[0] != '#'){
            std::istringstream iss(line);
            double x, x1;
            iss >> x >> x1;
            s.push_back(x);
            n.push_back(x1);
        }
    }
}

double linear_interpolate(double x1, double y1, double x2, double y2, double x){
    if (x1 > x2) {
        std::cout << "linear interpolation: points may be in the wrong order: x1 "<<  x1 << " x2 " << x2 << std::endl;
    }
    double a = ( y2 - y1 ) / ( x2 - x1 );
    double b = y2 - a * x2;
    return a * x + b;
}

std::vector<double> calculate_stable(const std::vector<double>& n, const std::vector<double>& s, const std::vector<double>& nx){
    std::vector<double> out;
    for (unsigned int i = 0; i < nx.size(); ++i){
        unsigned int x1 = left_index_min_zero_max_smallerlast(n, nx[i]);
        out.push_back(linear_interpolate(n[x1], s[x1], n[x1+1], s[x1+1], nx[i]));
    }
    return out;
}

std::vector<double> nstable(int Nsp, int& Nmulti){
    std::vector<double> n;
    std::vector<double> s;
    load_data(n, s);
    double maximum = *std::max_element(n.begin(), n.end());
    double step = maximum / double(Nsp);
    Nmulti = std::floor(step);
    std::vector<double> nx = arange(step, maximum, step);
    auto out = calculate_stable(n, s, nx);
    for (unsigned int i = 0; i < nx.size(); i++){
        std::cout << out[i] << " " << nx[i] << std::endl;
    }
    return out;
}
