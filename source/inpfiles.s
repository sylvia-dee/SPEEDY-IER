# Script to link fortran units to input files; $1 = resolution (t21, t30, ..)

SB=../data/bc/$1/clim
SC=../data/bc/$1/anom
SH=../hflux	

 ln -s $SB/orog_lsm_alb.${1}.grd         fort.20
 ln -s $SB/sst_8190clim.${1}.sea.grd     fort.21
 ln -s $SB/seaice_8190clim.${1}.sea.grd  fort.22
# ln -s $SB/skt_8190clim.${1}.land.grd    fort.23
# ln -s $SB/surfv_st_8190clim.t30.land.grd   fort.23
 ln -s $SB/ccsm4_skint_8190clim.t30.global.grd  fort.23	
 ln -s $SB/sndep_8190clim.${1}.land.grd  fort.24
 ln -s $SB/veget.${1}.land.grd           fort.25
 ln -s $SB/soilw_8190clim.${1}.land.grd  fort.26

# cp    $SC/sst_anom_1950_1999_mean8190.${1}.grd  fort.30
# cp    $SC/sst_anom_7990.${1}.grd	         fort.30
# cp    $SC/noaa_anom_1854_2002_mean8190.${1}.grd fort.30
# cp    $SC/hadisst_anom_1_1_1870_2002_mean8190.${1}.grd fort.30		
# cp    $SC/sst_anom_8190.${1}.sea.grd	         fort.30
# cp    $SC/elnino_anom_p1.${1}.grd	         fort.30
 cp    $SC/ccsm4_ts_anom_1_1_850_2005_mean8190.t31.grd fort.30

cp    $SC/noaa_anom_1854_2002_mean8190.${1}.grd fort.40	
	
cp    $SH/hflux_a20_a29_1949_2002_clim.grd  fort.31
#cp    $SH/hflux_zero.grd  fort.31	

ln -s $SB/ocean_model_sst.grd  fort.50
