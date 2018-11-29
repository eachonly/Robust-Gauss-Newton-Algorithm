# Compile and run RGN hymod demo using ifort

# Compile
ifort -c ../../SRC_RGN/constantsMod.f90
ifort -c ../../SRC_RGN/rgn.f90
ifort -c ../../SRC_DEMO/hymod/hymod.f90
ifort -c ../../SRC_DEMO/hymod/hydroDataMod.f90
ifort -c ../../SRC_DEMO/hymod/rgnMain_Hymod.f90

# Link
ifort constantsMod.o rgn.o hymod.o hydroDataMod.o rgnMain_Hymod.o -o hymod_demo.exe

echo "Build step done"
