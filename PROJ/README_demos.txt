The instructions below are for building rosenbrock demo
To use a different demo, replace "rosenbrock" with demo name, e.g. "hymod"


A. To build rosenbrock demo with ifort in Windows

Option 1: Use batch file rosenbrock_ifort.bat
Option 2: Use IDE, e.g. Visual Studio


B. To build rosenbrock demo with gfortran in Windows

Option 1: Use batch file rosenbrock_gfortran.bat
Option 2: Use IDE, e.g. Codeblocks


C. To build rosenbrock demo with ifort in Linux

Option 1: Use makefile
      make -f makefile.rosenbrock_ifort
      make -f makefile.rosenbrock_ifort clean # optional to remove intermediates
      ./rosenbrock_demo.exe

Option 2: Use bash file
      bash ./rosenbrock_ifort.sh
      ./rosenbrock_demo.exe

Note that the bash file can be generated automatically from the makefile
      make -f makefile.rosenbrock_ifort -n >rosenbrock_ifort.sh


D. To compile and run rosenbrock demo with gfortran in Linux

Same commands as for (C), but replace ifort with gfortran in filenames

