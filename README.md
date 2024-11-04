This is a solver for bending static analysis of a rectangular body with one face fixed and loads applied to its opposing faces.\
This program have two part, mesher for rectangular body (see below, excution example) and solver for static analysis.\
Static analysis solver can solve small deformation analysis (linear) and large deformation analysis (nonlinear).\
Newton-Raphson method as a non-linear analysis is implemented in this program. \
You can choose these analysis by setting a flag in an input file (./example/input.dat).\
You can also change Young modulus, Poisson ratio, density, and maximum iteration of Newton-Raphson method in the input file.
## install manual for linux
You need following environments.
1. make
2. cmake
3. git
4. bcc/gfortran
5. Lapack
6. MPI

**commnd example for install**
```
git clone git@github.com:HiranoAkihiro/FEM_large_deformation_analysis.git
./monolis_install.sh
```
## excution example
```
./mesher.sh 50.0 1.0 1.0 200 4 4
make
./run.sh
```
The part marked `./mesher.sh 50.0 1.0 1.0 200 4 4` should be entered as \
`./mesher.sh A B C D E F`

- A is length in x-direction.
- B is length in y-direction.
- C is length in z-direction.
- D is number of divisions in x-direction.
- E is number of divisions in y-direction.
- F is number of divisions in z-direction.
