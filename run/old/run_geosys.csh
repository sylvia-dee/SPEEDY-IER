#!/bin/csh 

# ${RES} = resolution (eg t21, t30)
# ${EXP} = experiment no. (eg 111)
# ${REST} = experiment no. for restart file ( 0 = no restart ) 

# Define directory names
# set -x

#if ( $# != 3 ) then
#    echo 'Usage: '$0' version resol. exp_no. restart_no' 1>&2
#    exit 1
#endif

set RES  = 't30'
set EXP  = 'dcn'
set REST = 0

set usr_src = 'nothing'


set UT=..
set SA=../source	
set CA=../tmp
set CB=../output/exp_${EXP}
set CC=../input
set CD=../output/exp_${REST}
set SAVDAT=../results	

mkdir ../output/exp_${EXP}	
		

# Edit input files if needed and store them with exp code

##echo "Do you want to modify the time-stepping parameters (y/n)?"
##read MSTEP
##if [ $MSTEP = 'y' ] ; then
##  vi $UT/ver40.input/cls_instep.h $SA/doc_instep.txt
##fi
##
##echo "Do you want to modify the dynamics parameters (y/n)?"
##read MSTEP
##if [ $MSTEP = 'y' ] ; then
##  vi $UT/ver40.input/cls_indyns.h  $SA/doc_indyns.txt
##fi
##
##echo "Do you want to modify the physics parameters (y/n)?"
##read MSTEP
##if [ $MSTEP = 'y' ] ; then
##  vi $UT/ver40.input/cls_inphys.h  $SA/doc_inphys.txt
##fi
##
##echo "Do you want to modify the input files (y/n)?"
##read MSTEP
##if [ $MSTEP = 'y' ] ; then
##  vi $UT/ver40.input/inpfiles.s
##fi


mkdir $UT/input/exp_${EXP}

echo "hor. resolution : " ${RES} >> $UT/input/exp_${EXP}/run_setup
echo "experiment no.  : " ${EXP} >> $UT/input/exp_${EXP}/run_setup
echo "restart exp. no.: " ${REST} >> $UT/input/exp_${EXP}/run_setup
	
# Copy files from basic version directory

pwd
echo "copying from $SA to $CA"
rm -f $CA/*

cp $SA/makefile $CA/
cp $SA/*.f      $CA/
cp $SA/*.h      $CA/
cp $SA/*.s      $CA/

cp $CA/par_horres_${RES}.h         $CA/atparam.h
#cp $CA/par_verres.h         $CA/atparam1.h
cp ../${usr_src}/par_verres.h      $CA/atparam1.h 

# Copy parameter and namelist files from user's .input directory

echo "ver40.input new files ..."
ls $UT/ver40.input

echo "copying parameter and namelist files from $UT/ver40.input "
cp $UT/ver40.input/cls_*.h     $CA/
cp $UT/ver40.input/inpfiles.s  $CA/
cp $UT/ver40.input/cls_*.h     $UT/input/exp_${EXP}
cp $UT/ver40.input/inpfiles.s  $UT/input/exp_${EXP}

# Copy modified model files from user's update directory

echo "update new files ..."
ls $UT/${usr_src}

echo "-----------------------------------------"
echo "copying modified model files from $UT/${usr_src}"
cp $UT/${usr_src}/* $CA/
cp $UT/$usr_src/* $UT/input/exp_${EXP}
	
# Set input files

cd $CA

# Set experiment no. and restart file (if needed)

echo ${REST} >  fort.2
echo ${EXP} >> fort.2

if ( ${REST} != 0 ) then
  echo "link restart file atgcm${REST}.rst to fort.3"
  ln -s $CD/atgcm${REST}.rst fort.3
endif

# Link input files

echo 'link input files to fortran units'

echo "-----------------------------------------"
pwd
bash inpfiles.s ${RES}
echo "-----------------------------------------"



ls -l fort.*

echo ' compiling at_gcm - calling make'

pwd
make imp.exe 
	


#
# create and execute a batch job to run the model
#

cat > run.job << EOF1
#QSUB -lM $5 -lT $6 -mb -me -r speedy -s /bin/csh -l mpp_p=1
 
#limit stacksize 150000

cd $CA
pwd

#
###---modification for coupling start		
#
 
mkdir ${SAVDAT}
/bin/rm ${SAVDAT}/fort.*	

#  cpl.org file is the fluxes for the first ocean day on the atmos model grid, if this file does not exist model fails	
cp ${SAVDAT}/cpl.org ${SAVDAT}/fort.105
	
#   Prepare data signal controlers for writer and reader
cat > ${SAVDAT}/fort.101 << --
0	
--

cp ${SAVDAT}/fort.101 ${SAVDAT}/fort.104
	
cat > ${SAVDAT}/fort.103 << --
1
--

cp ${SAVDAT}/fort.103 ${SAVDAT}/fort.106

#
# fort.200 is character string directy where the results are saved and both ocean and atmos
# write and read from
#
	
cat > fort.200 << --
'${SAVDAT}'	
--
		
#
# junk generate fort.102 which is "empty" (i.e., no sst ready to be read)
# ocean will update fort.102 when sst is ready
	
/bin/rm junk.f
cat > junk.f << --

        program main
        real*8 fnum
        fnum=102.0
        write(102)fnum
        stop
        end
--
ifort -convert big junk.f
a.out
mv fort.102 ${SAVDAT}/fort.102
/bin/rm a.out			

#		
###---modification for coupling end
#
				
setenv F_UFMTENDIAN big
#export F_UFMTENDIAN

echo 'the executable file...'
ls -l imp.exe


#
# RUN MODEL
#
	
set stacksize unlimited
set datasize unlimited

time ./imp.exe | tee out.lis

mv out.lis $CB/atgcm${EXP}.lis
mv fort.10 $CB/atgcm${EXP}.rst

mv at*${EXP}.ctl   $CB
mv at*${EXP}_*.grd $CB

mv day*${EXP}.ctl   $CB
mv day*${EXP}_*.grd $CB

cp fort.51 $CB	
	
cd $CB

chmod 644 at*${EXP}.*

#
###---modification for coupling start		
#
#  SAVE DATA
#
echo "saving results and restart data" 
cp ${SAVDAT}/fort.105 ${SAVDAT}/cpl
#
###---modification for coupling start		
#
EOF1

#qsub run.job

csh run.job 

exit
