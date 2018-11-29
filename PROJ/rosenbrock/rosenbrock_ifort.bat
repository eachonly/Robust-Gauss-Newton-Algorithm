@echo off

:: Compile and run RGN rosenbrock demo using ifort
echo Compile rosenbrock demo with ifort

:: Ensure correct PATH for ifort
set "ifortvars_path=C:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\"
echo Assuming path to ifort is %ifortvars_path%
set /p input=Is this the correct installation of Intel Fortran? (Y/N) 

if /i "%input%"=="Y" (
:: To build 32 bit executable replace "intel64" with "ia32"
  call "%ifortvars_path%ifortvars.bat" intel64
) else (
  echo Please modify .bat to have correct path to ifortvars and retry
  pause
  exit /b
)

:: Compile
ifort -c ..\..\SRC_RGN\constantsMod.f90
ifort -c ..\..\SRC_RGN\rgn.f90
ifort -c ..\..\SRC_DEMO\rosenbrock\rgnMain_Rosenbrock.f90

:: Link
ifort -o rosenbrock_demo.exe *.obj

echo Build step done

:: Run
.\rosenbrock_demo.exe

pause
