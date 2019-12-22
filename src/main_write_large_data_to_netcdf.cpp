#include "logger.h"
#include "yaml-cpp/yaml.h"
#include "setupcolumnmodelyaml.h" 
#include <vector>
#include "superparticle.h"


auto configdata = YAML::Load( R"FOO(
t_max: 3000
dt: 0.1
grid:
    toa: 3000.
    gridlength: 1.
initial_state:
    ALR: 0.004
    Theta0: 297.2
    p0: 100000.
    cloud_base: 500.
    cloud_roof: 500.
    w: 2.
)FOO");

int main(){
    auto logger = createLogger("netcdf", "out_write_large_data_to_netcdf.nc");
    auto grid = createGrid(configdata["grid"]);
    auto state = createState(*grid, configdata["initial_state"]);
    logger->initialize(state, 0.1);
    std::vector<Superparticle> sps(300);
    for (int i=0; i<30000; ++i)
    {
        logger->log(state, sps);
    }
}
