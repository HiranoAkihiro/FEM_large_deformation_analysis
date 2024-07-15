#!/bin/bash
cd submodule
git submodule add git@github.com:nqomorita/monolis.git
cd monolis
./install.sh
make
cd ../../