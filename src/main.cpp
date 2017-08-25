#include <iostream>
#include "columnmodel.h"
//#include "setupcolumnmodelyaml.h"
#include "setupcolumnmodeldummy.h"
#include "logger.h"
#include "time_stamp.h"

int main(int argc, char** argv) {
//    if (argc < 2) {
//        std::cerr << "Need Yaml Filename" << std::endl;
//        return -1;
//    }
//    try {
//        const YAML::Node config = YAML::LoadFile(argv[1]);
        auto columnmodel = createColumnModel();
        std::string file_name = "dummy.nc";
        //std::string file_name = time_stamp();
        std::shared_ptr<Logger> logger = createLogger("netcdf", file_name);
        columnmodel.run(logger);
//    }
//    catch (YAML::Exception e) {
//        std::cerr << "Error while parsing yaml file" << std::endl;
//        std::cerr << "Line: " << e.mark.line << " Col: " << e.mark.column
//                  << std::endl;
//        std::cerr << e.what() << std::endl;
//        throw e;
//    } catch (std::invalid_argument e) {
//        std::cerr << e.what() << std::endl;
//    }
}
