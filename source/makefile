#
# This makefile for running on catfish.colorado.edu
# David Noone <dcn@colorado.edu> - Fri Feb 18 19:21:08 MST 2011
# #

COMPILE=ifort 
COMLIB1=

#
# Compiler flags listed with most optomize/fastest at top, to 
# most diagnostic/slowest at the bottom.
#
#COMOTT1=  -ipo -O3 -no-prec-div -static
#COMOTT1=  -O3 -no-prec-div -static
COMOTT1=  -O1 -fpe0 -fp-stack-check -traceback -g 
#COMOTT1=  -O1 -CB -CU -fpe0 -fp-stack-check -traceback -g
#COMOTT1 = -check all -debug all -fpe0 -fp-stack-check -traceback -g 

 COMCONV= -convert big 
# COMCONV=

FILES= \
       addflx.o \
       atm2land.o \
       atm2sea.o \
       get_from_coupler.o \
       dyn_geop.o \
       dyn_grtend.o \
       dyn_implic.o \
       dyn_sptend.o \
       dyn_step.o \
       dyn_stloop.o \
       ini_coupler.o \
       ini_impint.o \
       ini_indyns.o \
       ini_inforc.o \
       ini_iniall.o \
       ini_inphys.o \
       ini_iniso.o \
       ini_inirdf.o \
       ini_invars.o \
       ini_stepone.o \
       land2atm.o \
       land_init.o \
       phy_convmf.o \
       phy_culnd.o  \
       phy_fordate.o \
       phy_lscond.o \
       phy_pevap.o \
       phy_phypar.o \
       phy_radiat.o \
       phy_shtorh.o \
       phy_suflux.o \
       phy_vdifsc.o \
       ppo_diagns.o \
       ppo_restart.o \
       ppo_setctl.o \
       ppo_setctl_daily.o \
       ppo_setgrd.o \
       ppo_tminc.o \
       ppo_tminc_daily.o \
       ppo_tmout.o \
       ppo_tmout_daily.o \
       iso_alpice.o \
       iso_alpliq.o \
       iso_alpkoc.o \
       iso_alpkss.o \
       iso_alpkex.o \
       iso_conden.o \
       iso_rescalr.o \
       iso_delta3d.o \
       iso_delta2d.o \
       iso_chkq.o \
       sea2atm.o \
       sea_init.o \
       send_to_coupler.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o \
       readerA.o writerA.o

.SUFFIXES: .f .o

.f.o:
	$(COMPILE)  $(COMOTT1) $(COMCONV) -c  $<

speedy.exec: $(FILES)  at_gcm.o
	$(COMPILE)  $(COMOTT1) $(COMCONV)  $(FILES) at_gcm.o -o speedy.exec $(COMLIB1) 


at_gcm.o               : at_gcm.f
addflx.o               : addflx.f
atm2land.o             : atm2land.f
atm2sea.o              : atm2sea.f
get_from_coupler.o     : get_from_coupler.f
dyn_geop.o             : dyn_geop.f
dyn_grtend.o           : dyn_grtend.f
dyn_implic.o           : dyn_implic.f
dyn_sptend.o           : dyn_sptend.f
dyn_step.o             : dyn_step.f
dyn_stloop.o           : dyn_stloop.f
ini_coupler.o          : ini_coupler.f
ini_impint.o           : ini_impint.f
ini_indyns.o           : ini_indyns.f
ini_inforc.o           : ini_inforc.f
ini_iniall.o           : ini_iniall.f
ini_inphys.o           : ini_inphys.f
ini_inirdf.o           : ini_inirdf.f
ini_invars.o           : ini_invars.f
ini_stepone.o          : ini_stepone.f
land2atm.o             : land2atm.f
land_init.o            : land_init.f
phy_convmf.o           : phy_convmf.f
phy_culnd.o            : phy_culnd.f
phy_lscond.o           : phy_lscond.f
phy_pevap.o            : phy_pevap.f
phy_fordate.o          : phy_fordate.f
phy_phypar.o           : phy_phypar.f
phy_radiat.o           : phy_radiat.f
phy_shtorh.o           : phy_shtorh.f
phy_suflux.o           : phy_suflux.f
phy_vdifsc.o           : phy_vdifsc.f
ppo_diagns.o           : ppo_diagns.f
ppo_restart.o          : ppo_restart.f
ppo_setctl.o           : ppo_setctl.f
ppo_setctl_daily.o     : ppo_setctl_daily.f
ppo_setgrd.o           : ppo_setgrd.f
ppo_tminc.o            : ppo_tminc.f
ppo_tminc_daily.o      : ppo_tminc_daily.f
ppo_tmout.o            : ppo_tmout.f
ppo_tmout_daily.o      : ppo_tmout_daily.f
iso_alpice.o           : iso_alpice.f
iso_alpliq.o           : iso_alpliq.f
iso_alpkoc.o           : iso_alpkoc.f
iso_alpkss.o           : iso_alpkss.f
iso_alpkex.o           : iso_alpkex.f
iso_conden.o           : iso_conden.f
iso_rescalr.o          : iso_rescalr.f
iso_delta3d.o          : iso_delta3d.f
iso_delta2d.o          : iso_delta2d.f
iso_chkq.o             : iso_chkq.f
sea2atm.o              : sea2atm.f
sea_init.o             : sea_init.f 
send_to_coupler.o      : send_to_coupler.f
spe_matinv.o           : spe_matinv.f
spe_spectral.o         : spe_spectral.f
spe_subfft_fftpack.o   : spe_subfft_fftpack.f 
writerA.o              : writerA.f
readerA.o              : readerA.f
