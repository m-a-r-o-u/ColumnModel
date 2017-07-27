#!/bin/bash

export CMAKE_PREFIX_PATH=$HOME/Software/yaml-cpp/install/';'
export CMAKE_PREFIX_PATH=$HOME/Software/rrtm/';'${CMAKE_PREFIX_PATH}
echo ${CMAKE_PREFIX_PATH}
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
