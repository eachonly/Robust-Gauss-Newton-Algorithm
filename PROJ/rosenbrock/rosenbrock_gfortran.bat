@echo off

:: Compile and run RGN rosenbrock demo using gfortran
echo Compile rosenbrock demo with gfortran

:: The version of executable file (32 bit or 64 bit) depends on the gfortran version

:: Ensure correct PATH for gfortran
set /p input=Is gfortran on the path environment variable? (Y/N) 

if /i "%input%"=="N" (
  echo Assuming PATH=C:\cygwin64\bin\ (otherwise, update .bat file)
  PATH=%PATH%;C:\cygwin64\bin\
)

:: Compile
gfortran -c -ffree-line-length-none ..\..\SRC_RGN\constantsMod.f90 ..\..\SRC_RGN\rgn.f90 ..\..\SRC_DEMO\rosenbrock\rgnMain_Rosenbrock.f90

:: Link
gfortran -o rosenbrock_demo.exe *.o

echo Build step done

:: Run
.\rosenbrock_demo.exe

pause

