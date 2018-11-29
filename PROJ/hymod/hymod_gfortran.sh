# Compile and run RGN hymod demo using gfortran

# Compile
gfortran -c -ffree-line-length-none ../../SRC_RGN/constantsMod.f90
gfortran -c -ffree-line-length-none ../../SRC_RGN/rgn.f90
gfortran -c -ffree-line-length-none ../../SRC_DEMO/hymod/hymod.f90
gfortran -c -ffree-line-length-none ../../SRC_DEMO/hymod/hydroDataMod.f90
gfortran -c -ffree-line-length-none ../../SRC_DEMO/hymod/rgnMain_Hymod.f90

# Link
gfortran constantsMod.o rgn.o hymod.o hydroDataMod.o rgnMain_Hymod.o -o hymod_demo.exe

echo "Build step done"
