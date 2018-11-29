# Compile and run RGN rosenbrock demo using ifort

# Compile
ifort -c ../../SRC_RGN/constantsMod.f90
ifort -c ../../SRC_RGN/rgn.f90
ifort -c ../../SRC_DEMO/rosenbrock/rgnMain_Rosenbrock.f90

# Link
ifort constantsMod.o rgn.o rgnMain_Rosenbrock.o -o rosenbrock_demo.exe

echo "Build step done"
