#!/bin/bash

gfortran example/mesher/mesher.f90 -o example/mesher/mesh
./example/mesher/mesh $1 $2 $3 $4 $5 $6