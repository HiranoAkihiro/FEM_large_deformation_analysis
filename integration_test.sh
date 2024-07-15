#!/bin/bash
cd directory_for_test
make clean
make
cd ../
./bin/integration_test