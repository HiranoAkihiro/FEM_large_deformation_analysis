### install manual for linux
You need following environments.
1. make
2. cmake
3. git
4. bcc/gfortran
5. Lapack
6. MPI

**commnd example**
```
git clone git@github.com:HiranoAkihiro/FEM_large_deformation_analysis.git
./monolis_install.sh
```
### excution example
```
./mesher.sh 50.0 1.0 1.0 200 4 4
make
./run.sh
```
After `./mesher.sh`, enter (length in x-direction) (length in y-direction) (length in z-direction) (number of divisions in x-direction) (number of divisions in y-direction) (number of divisions in z-direction) in that order.
