#!/bin/csh -f

#my job to launch
set myJob=Slurmjob_Bsp2D.sh

#my directory for compiling
set mydirComp=Bspline-2D

#my directory for running jobs 
set mydirRun=$USERWORK/Bspline2D

#copy the slurm job to mydid
cp "${myJob}" "${mydirRun}"

echo --This script runs the codes on:
echo ----"${mydirRun}"----

set exec = bspline_test-2d.exe

cd ../${mydirComp}
#check if the executable is generated before lauching the Slurm scrip
if ( -e ${exec} ) then
   echo ---The executable is generated--
else
  echo --The file -"${exec}"- does not exist-- I stop
  exit
endif

#copy the executable to mydir
cp "${exec}" "${mydirRun}"

cd "${mydirRun}"

#launch the Slurmjob
sbatch ${myJob}
