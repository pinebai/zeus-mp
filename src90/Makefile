#=======================================================================
#
#
#  Linux MAKEFILE FOR ZEUS-MP
#
#  Contributed by Asif ud-Doula
#
#  NOTE: This Makefile does not include links to FFTw 2.1.5 system
#        libraries.  If the FFT solver is to be used, the following
#        changes to this Makefile must be made:
#        (1) appropriate paths for linking to the FFTw 2.1.5 library
#            must be added.
#        (2) the string "-DFFT" must be added to the ZMP_CPP definition
#            below (line 40).
#        (3) fftw_ps.o and fftwplan.o must be added to the COBJ definition
#            (line 67).
#
#-------------------  macro definitions  -------------------------------
#
all: directories compile
ZMP_EXE   = ../exe90/
ZMP_MOD   = ./
ZMP_FC    = mpifort -ffixed-line-length-0
MKDIR_P = mkdir -p
#-----------------------------------------------------------------------
# Uncomment the next line to compile with agressive optimization
#ZMP_OPTS  = -c -O3
ZMP_OPTS  = -c -g
#-----------------------------------------------------------------------
ZMP_CC    = mpicc
ZMP_LDR   = mpifort -ffixed-line-length-0
#-----------------------------------------------------------------------
ZMP_CPP   = cpp -DPROBLEM=blast2d -DPROBRES=empty -DMPI_USED -DMPI_FLOAT=MPI_DOUBLE_PRECISION -DMPI_2FLOAT=MPI_2DOUBLE_PRECISION
#--------------------  implicit rules  ---------------------------------
#
.SUFFIXES: .F .for .c .o

.for.o:
	${ZMP_CPP} -P -traditional $*.F | tr -s "" "\012" > $*.for
	${ZMP_FC} -I ${ZMP_MOD} ${ZMP_OPTS} $*.for

.F.o:
	${ZMP_CPP} -P -traditional $*.F | tr -s "" "\012" > $*.for
	${ZMP_FC} -I ${ZMP_MOD} ${ZMP_OPTS} $*.for

.c.o:
	${ZMP_CC} -c $<

#
#-------------------  Fortran object files  ----------------------------
#
FOBJ = mod_files.o advx1.o advx2.o advx3.o alloc_arrays.o avisc.o avisc_d.o binout.o bndyflgs.o bval3d.o bvalabuns.o bvalemf.o bvalt.o cgsolve.o crdiffuse.o crdiffuse_d.o crsource.o clocks.o configure.o cooling.o ct.o ct_1D.o ct_2D.o dataio.o difco1.o difco2.o difco3.o diverg.o dprd.o empty.o eos_d.o eos.o fftwgrav.o findno.o floors.o forces.o forces_d.o ggen.o gpbv.o gravity.o gravmg.o grav2dcg.o grav3dcg.o grdv.o hdfall.o hsmoc.o hsmoc_1D.o hsmoc_2D.o grey_fld.o intchk.o jacob.o lorentz.o lorentz_1D.o lorentz_2D.o lorentz_d.o mapin.o mapout.o matprop.o mnmx.o momx1.o momx2.o momx3.o msave.o msendrec.o mstart.o mul_sym.o offdiag.o newdt.o newgrid.o newvg.o newx1.o newx2.o newx3.o nudt.o opac_d.o opacity.o options.o output2D.o pdv.o pdv_d.o printd.o qupdate.o qupdate_d.o qupdate_iso.o qupdate_iso_d.o rad_imp_dt.o rad_solve.o radialplot.o restart.o sendrec.o setup.o srcstep.o strtoi.o sym_diag.o sym_prcn.o transprt.o transprt_1D.o transprt_2D.o tranx1.o tranx2.o tranx3.o updtmtrxbnd.o updtvecbnd.o zeusmp.o gcollapse.o textdmp.o  linpck.o sod.o  sedov.o stream.o rshock.o mhdshktube.o rotor.o mhdvortex.o noh.o torus.o cluster.o bondi.o thermal.o blast2d.o 
#
#-------------------  C object files  ----------------------------------
#
COBJ = bcdflt.o checkin.o
#
#---------------------  targets  ---------------------------------------
#
help:
	@echo Type "'make help'             to see this menu"
	@echo Type "'make compile           to create executable
	@echo Type "'make clean'            to remove '*.for' and '*.o' files"
	@echo Type "'make cleaner'          to remove '*.for' '*.o' "
	@echo "                             ${ZMP_EXE}zmp_inp ${ZMP_EXE}zeusmp.x"
	@echo Type "'make listing'          to create listing of code"
#
#-----------------------------------------------------------------------
#
NPFOBJ = setup.for restart.for textdmp.for forces.for hdfall.for opacity.for matprop.for bval3d.for difco1.for difco2.for difco3.for
NPOOBJ = setup.o restart.o textdmp.o forces.o hdfall.o opacity.o matprop.o bval3d.o difco1.o difco2.o difco3.o

newprob:
	\rm -f ${NPFOBJ} ${NPOOBJ}
#
#-----------------------------------------------------------------------
#
clean:
	echo "\rm -f *.o *.for *.mod"
	\rm -f *.o *.for *.mod
#
#-----------------------------------------------------------------------
#
cleaner:
	echo "\rm -f *.o *.for *.mod  ${ZMP_EXE}zmp_inp ${ZMP_EXE}zeusmp.x"
	\rm -f *.o *.for  ${ZMP_EXE}zmp_inp ${ZMP_EXE}zeusmp.x
#
#-----------------------------------------------------------------------
#
listing:
	echo "cat Makefile  *.h *.F *.c > listing"
	cat Makefile  *.h *.F *.c > listing
#
#-----------------------------------------------------------------------
#
directories:${ZMP_EXE}
${ZMP_EXE}:
	${MKDIR_P} ${ZMP_EXE}

compile:${ZMP_EXE}zeusmp.x

${ZMP_EXE}zeusmp.x: ${FOBJ} ${COBJ}
	${ZMP_LDR} -o ${ZMP_EXE}zeusmp.x ${FOBJ} ${COBJ}

#
#----------------------  dependencies  ---------------------------------
#
mod_files.o: mod_files.F 

advx1.o: advx1.F mod_files.o  

advx2.o: advx2.F mod_files.o  

advx3.o: advx3.F mod_files.o  

alloc_arrays.o: alloc_arrays.F mod_files.o

avisc.o: avisc.F mod_files.o  

avisc_d.o: avisc_d.F mod_files.o  

binout.o: binout.F mod_files.o

bndyflgs.o: bndyflgs.F mod_files.o  

bval3d.o: bval3d.F mod_files.o  

bvalabuns.o: bvalabuns.F mod_files.o  

bvalemf.o: bvalemf.F mod_files.o  

bvalt.o: bvalt.F mod_files.o  

cgsolve.o: cgsolve.F mod_files.o

clocks.o: clocks.F mod_files.o  

configure.o: configure.F mod_files.o  

cooling.o: cooling.F mod_files.o

crdiffuse.o: crdiffuse.F mod_files.o

crdiffuse_d.o: crdiffuse_d.F mod_files.o

crsource.o: crsource.F mod_files.o

ct.o: ct.F mod_files.o  

ct_1D.o: ct_1D.F mod_files.o  

ct_2D.o: ct_2D.F mod_files.o  

dataio.o: dataio.F mod_files.o  

difco1.o: difco1.F mod_files.o

difco2.o: difco2.F mod_files.o

difco3.o: difco3.F mod_files.o

diverg.o: diverg.F mod_files.o  

dprd.o: dprd.F mod_files.o

empty.o: empty.F  

eos_d.o: eos_d.F mod_files.o

eos.o: eos.F mod_files.o

findno.o: findno.F  

floors.o: floors.F mod_files.o

forces.o: forces.F mod_files.o  

forces_d.o: forces_d.F mod_files.o  

fourn.o: fourn.F  

ggen.o: ggen.F mod_files.o  

gravity.o: gravity.F mod_files.o  

gravmg.o: gravmg.F mod_files.o  

grav2dcg.o: grav2dcg.F mod_files.o  

grav3dcg.o: grav3dcg.F mod_files.o  

grdv.o: grdv.F mod_files.o

gpbv.o: gpbv.F mod_files.o

grdv.o: grdv.F mod_files.o  

hdfall.o: hdfall.F mod_files.o  

hsmoc.o: hsmoc.F mod_files.o  

hsmoc_1D.o: hsmoc_1D.F mod_files.o  

hsmoc_2D.o: hsmoc_2D.F mod_files.o  

grey_fld.o: grey_fld.F mod_files.o

intchk.o: intchk.F mod_files.o  

jacob.o: jacob.F mod_files.o

lorentz.o: lorentz.F mod_files.o  

lorentz_1D.o: lorentz_1D.F mod_files.o  

lorentz_2D.o: lorentz_2D.F mod_files.o  

lorentz_d.o: lorentz_d.F mod_files.o  

mapin.o: mapin.F mod_files.o  

mapout.o: mapout.F mod_files.o  

matprop.o: matprop.F mod_files.o

msendrec.o: msendrec.F mod_files.o

mnmx.o: mnmx.F mod_files.o  

momx1.o: momx1.F mod_files.o  

momx2.o: momx2.F mod_files.o  

momx3.o: momx3.F mod_files.o  

msave.o: msave.F mod_files.o  

mstart.o: mstart.F mod_files.o  

mul_sym.o: mul_sym.F mod_files.o

newdt.o: newdt.F mod_files.o  

newgrid.o: newgrid.F  mod_files.o

newvg.o: newvg.F mod_files.o

newx1.o: newx1.F mod_files.o

newx2.o: newx2.F mod_files.o

newx3.o: newx3.F mod_files.o

nudt.o: nudt.F mod_files.o  

offdiag.o: offdiag.F mod_files.o

opac_d.o: opac_d.F mod_files.o

opacity.o: opacity.F mod_files.o

options.o: options.F mod_files.o 

output2D.o: output2D.F mod_files.o

pdv.o: pdv.F mod_files.o  

pdv_d.o: pdv_d.F mod_files.o  

printd.o: printd.F mod_files.o

qupdate.o: qupdate.F mod_files.o

qupdate_d.o: qupdate_d.F mod_files.o

qupdate_iso.o: qupdate_iso.F mod_files.o

qupdate_iso_d.o: qupdate_iso_d.F mod_files.o

rad_imp_dt.o: rad_imp_dt.F mod_files.o

rad_solve.o: rad_solve.F mod_files.o

radialplot.o: radialplot.F mod_files.o

restart.o: restart.F mod_files.o  

sendrec.o: sendrec.F mod_files.o

setup.o: setup.F mod_files.o  

srcstep.o: srcstep.F mod_files.o  

strtoi.o: strtoi.F  

sym_diag.o: sym_diag.F mod_files.o

sym_prcn.o: sym_prcn.F mod_files.o

tmul_sym.o: tmul_sym.F mod_files.o

transprt.o: transprt.F mod_files.o  

transprt_1D.o: transprt_1D.F mod_files.o  

transprt_2D.o: transprt_2D.F mod_files.o  

tranx1.o: tranx1.F mod_files.o  

tranx2.o: tranx2.F mod_files.o  

tranx3.o: tranx3.F mod_files.o  

updtmtrxbnd.o: updtmtrxbnd.F mod_files.o  

updtvecbnd.o: updtvecbnd.F mod_files.o  

zeusmp.o: zeusmp.F mod_files.o  

textdmp.o: textdmp.F mod_files.o  

linpck.o: linpck.F mod_files.o  

bcdflt.o: bcdflt.c  

checkin.o: checkin.c  

rotor.o: rotor.F mod_files.o  

mhdvortex.o: mhdvortex.F mod_files.o  

mhdshktube.o: mhdshktube.F mod_files.o  

blast2d.o: blast2d.F mod_files.o  

gcollapse.o: gcollapse.F mod_files.o  

sedov.o: sedov.F mod_files.o  

sod.o: sod.F mod_files.o  

stream.o: stream.F mod_files.o

torus.o: torus.F mod_files.o

cluster.o: cluster.F mod_files.o

thermal.o: thermal.F mod_files.o

blast.o: blast.F mod_files.o

#shock_tube.o: shock_tube.F mod_files.o

bondi.o: bondi.F mod_files.o

rshock.o: rshock.F mod_files.o

noh.o: noh.F mod_files.o

fftw_ps.o: fftw_ps.c

fftwplan.o: fftwplan.c
