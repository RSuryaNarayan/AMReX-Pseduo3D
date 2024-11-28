# AMReX-Pseduo3D
Takes in a 2-D AMReX plotfile and converts it to a 3D plotfile which is one cell thick along the third direction. 

Compile with [amrex](https://github.com/AMReX-Codes/amrex) and run. Tested inside  [```PeleAnalysis/Src```](https://github.com/AMReX-Combustion/PeleAnalysis/tree/master/Src) by using:
```
make -j<cores> EBASE=Pseudo3D NEEDS_f90_SRC=FALSE
```
Runtime inputs: 

```
./Pseudo3D3d.gnu.MPI.ex infile=<pltfilename> is_per=<i j 1>
```
For the second argument, we specifify periodicity of boundaries along the 2 directions. The third direction is kept periodic after extrusion. 

# Extrude3D 

Script meant to take in a 2-D AMReX plotfile and extrude it uniformly along the third direction with a specified number of ```n_cells_z``` along the third direction. Would be a useful capability to speed up three-dimensional simulations as it gives a better initial condition that having to evolve things from the start in most cases. Currently doesn't work and needs fixing! 

