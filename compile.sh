#!/bin/bash

mkdir -p build && cd build && rm -rf ./*
cmake -DBUILD_TESTS=True -DBUILD_DOCS=True -DMPI=True -DVERBOSE=True -DCMAKE_BUILD_TYPE=Release ..
make
cd -
