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
The part marked `./mesher.sh 50.0 1.0 1.0 200 4 4` should be entered as \
`./mesher.sh A B C D E F`

- A is length in x-direction.
- B is length in y-direction.
- C is length in z-direction.
- D is number of divisions in x-direction.
- E is number of divisions in y-direction.
- F is number of divisions in z-direction.
