#include <iostream>
#include "columnmodel.h"
#include "setupcolumnmodelyaml.h"
//#include "setupcolumnmodeldummy.h"
#include <yaml-cpp/yaml.h>
#include "logger.h"
#include "time_stamp.h"

int main(int argc, char** argv) {
    try {
        YAML::Node config;
        if (argc == 2) {
            std::ifstream input(argv[1]);
            config = YAML::Load(input);
        } else {
            config = YAML::Load(std::cin);
        }
        std::random_device rd;
        std::mt19937_64 gen(rd());
        auto columnmodel = createColumnModel(gen, config["model"]);
        std::string file_name = "dummy.nc";
        std::shared_ptr<Logger> logger = createLogger(config["logger"]);
        columnmodel.run(logger);
    } catch (YAML::Exception e) {
        std::cerr << "Error while parsing yaml file" << std::endl;
        std::cerr << "Line: " << e.mark.line << " Col: " << e.mark.column
                  << std::endl;
        std::cerr << e.what() << std::endl;
        throw e;
    } catch (std::invalid_argument e) {
        std::cerr << e.what() << std::endl;
    }
}
