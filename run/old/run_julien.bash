#!/bin/bash

export LD_LIBRARY_PATH=/usr/geosys/x86_64/pgi/linux86-64/6.1/lib/
# $1 = resolution (eg t21, t30)
# $2 = experiment name (THREE CHARACTERS MAX) 
# $3 = experiment name for restart file ( 0 = no restart ) 

# Define directory names
# set -x

if [ $# -ne 3 ] ; then

    echo 'Usage: '$0' version resol. exp_no. restart_no' 1>&2
    exit 1

fi

UT=..
SA=../source	
CA=../tmp
mkdir ../output/exp_$2	
CB=../output/exp_$2
CC=../input
CD=../output/exp_$3

# Edit input files if needed and store them with exp code

echo "Do you want to modify the time-stepping parameters (y/n)?"
read MSTEP
if [ $MSTEP = 'y' ] ; then
  nedit $UT/ver40.input/cls_instep.h $SA/doc_instep.txt
fi

echo "Do you want to modify the dynamics parameters (y/n)?"
read MSTEP
if [ $MSTEP = 'y' ] ; then
  nedit $UT/ver40.input/cls_indyns.h  $SA/doc_indyns.txt
fi

echo "Do you want to modify the physics parameters (y/n)?"
read MSTEP
if [ $MSTEP = 'y' ] ; then
  nedit $UT/ver40.input/cls_inphys.h  $SA/doc_inphys.txt
fi

echo "Do you want to modify the input files (y/n)?"
read MSTEP
if [ $MSTEP = 'y' ] ; then
  nedit $UT/ver40.input/inpfiles.s
fi


mkdir $UT/input/exp_$2

echo "hor. resolution : " $1 >> $UT/input/exp_$2/run_setup
echo "experiment no.  : " $2 >> $UT/input/exp_$2/run_setup
echo "restart exp. no.: " $3 >> $UT/input/exp_$2/run_setup
	
# Copy files from basic version directory

echo "copying from $SA to $CA"
rm -f $CA/*

cp $SA/makefile $CA/
cp $SA/*.f      $CA/
cp $SA/*.h      $CA/
cp $SA/*.s      $CA/

mv $CA/par_horres_$1.h   $CA/atparam.h
mv $CA/par_verres.h      $CA/atparam1.h 

# Copy parameter and namelist files from user's .input directory

echo "ver40.input new files ..."
ls $UT/ver40.input

echo "copying parameter and namelist files from $UT/ver40.input "
cp $UT/ver40.input/cls_*.h     $CA/
cp $UT/ver40.input/inpfiles.s  $CA/
cp $UT/ver40.input/cls_*.h     $UT/input/exp_$2
cp $UT/ver40.input/inpfiles.s  $UT/input/exp_$2

# Copy modified model files from user's update directory

echo "update new files ..."
ls $UT/update
cp $SA/makefile $UT/update/

echo "copying modified model files from $UT/update"
cp $UT/update/* $CA/
cp $UT/update/* $UT/input/exp_$2
	
# Set input files

cd $CA

# Set experiment no. and restart file (if needed)

echo $3 >  fort.2
echo $2 >> fort.2

if [ $3 != 0 ] ; then
  echo "link restart file atgcm$3.rst to fort.3"
  ln -s $CD/atgcm$3.rst fort.3
fi 

# Link input files

echo 'link input files to fortran units'

bash inpfiles.s $1

ls -l fort.*

echo ' compiling isotope-enabled SPEEDY-IER - calling make'
pwd
make speedy.exec
	

#
# create and execute a batch job to run the model
# ==================================================
export PBS_0_WORKDIR=`pwd`
echo ' Assembling PBS script'
cat > run_speedy-ier.pbs << EOF1
#PBS -N SPEEDY
#PBS -q largemem
#PBS -o SPEEDY.out
#PBS -e SPEEDY.err
#PBS -S /bin/bash

export LD_LIBRARY_PATH=/opt/intel/lib/fce
export PBS_0_WORKDIR=$PBS_0_WORKDIR
export PATH=./:/usr/geosys/x86_64/intel/fce/9.0/bin:$PATH

cd $PBS_0_WORKDIR

set -x
 
#limit stacksize 150000

pwd

echo 'the executable file...'
ls -l speedy.exec

#
# RUN MODEL
#
	
ulimit -s unlimited
ulimit -d unlimited	

time ./speedy.exec > out.lis

mv out.lis $CB/atgcm$2.lis
mv fort.10 $CB/atgcm$2.rst

mv at*$2.ctl   $CB
mv at*$2_*.grd $CB

mv day*$2.ctl   $CB
mv day*$2_*.grd $CB
	
cd $CB

chmod 644 at*$2.*

EOF1

#qsub run_speedy-ier.pbs

exit
