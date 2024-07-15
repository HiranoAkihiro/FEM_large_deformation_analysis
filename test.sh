#!/bin/bash
cd directory_for_test/test
make clean
make
cd ../../
./bin/test