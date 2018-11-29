# Compile and run RGN rosenbrock demo using gfortran

# Compile
gfortran -c -ffree-line-length-none ../../SRC_RGN/constantsMod.f90
gfortran -c -ffree-line-length-none ../../SRC_RGN/rgn.f90
gfortran -c -ffree-line-length-none ../../SRC_DEMO/rosenbrock/rgnMain_Rosenbrock.f90

# Link
gfortran constantsMod.o rgn.o rgnMain_Rosenbrock.o -o rosenbrock_demo.exe

echo "Build step done"
