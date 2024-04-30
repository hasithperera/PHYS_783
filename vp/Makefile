#################################################################### OVERVIEW
#  Makefile for VP: 1D1V Vlasov-Poisson Nonlinear Solver
#
#   VERSION 3.0
#
#  Copyright, 2019
#  Gregory Howes, Kristopher Klein
#
#  Version notes:
#  KGK: Update to allow for arbitrary number of drifting Maxwellian distributions.
#  LAST UPDATE:  27 MAR 2019
###############################################################################

#COMPILERS AND FLAGS
#===========================================================
# INTEL: For Intel Fortran compiler ifort
#COMP=   ifort
#FLAGS=  -O4 -r8 -double-size 128

# OR DEBUG for ifort
#FLAGS=  -g -r8 -double-size 128

#===========================================================
# GFORTRAN:  For gfortran
COMP= gfortran
FLAGS= -O3 -DDOUBLE -ffree-form -ffree-line-length-512

# OR DEBUG for gfortran
#FLAGS= -g -DDOUBLE -ffree-form -ffree-line-length-512

#===========================================================
# IBM XL FORTRAN:  For IBM XL Fortran
#COMP=   xlf90
#FLAGS= -O4 -qautodbl=dbl -qsuffix=f=f90
#FLAGS=  -qautodbl=dbl -qsuffix=f=f90 -Iinclude

LIBS=	

VPMOD= 	vp_data.o vp_disp_funcs.o vp_funcs.o

VPX=  	vp.o

TESTVPX= test_vp.o

VPATH= 	src:include:/usr/include/

SH =    `which sh`

###############################################################################
all: clean vp neat

vp: $(VPMOD) $(VPX)
	$(COMP) -o vp.e $(FLAGS) $(LIBS) $(VPMOD) $(VPX)

test: 	 $(VPMOD) $(TESTVPX)
	$(COMP) -o test_vp.e $(FLAGS) $(LIBS) $(VPMOD) $(TESTVPX)

###############################################################################

###############################################################################

docs:	type bib type docs.ps neat_docs

pdf:	vp_docs.ps
	ps2pdf vp_docs.ps

type:	vp_docs.tex 
	${SH} -c "latex vp_docs.tex"

docs.ps:	vp_docs.dvi
	dvips -t letter -o  vp_docs.ps  vp_docs.dvi
	ps2pdf vp_docs.ps

bib:    
	echo "Attempting to make the bibliography.....";\
	${SH} -c "bibtex  vp_docs; latex  vp_docs; latex vp_docs"


neat:	
	mv *.mod include/
	mv *.o include/

neat_docs:
	mv vp_docs.aux docs/
	mv vp_docs.bbl docs/
	mv vp_docs.blg docs/
	mv vp_docs.dvi docs/
	mv vp_docs.log docs/
	mv vp_docs.out docs/
	mv vp_docs.toc docs/


clean:
	rm -f include/*.o include/*.mod
	rm -f vp*.e
	rm -f *.o

tar: 
	tar -cvf  vp3_`date +'%y%m%d'`.tar Makefile src/*.f90 include docs/*tex docs/*ps vp_docs.tex input *.in


#########Rules
%.o : %.f90
	$(COMP) -c $(FLAGS) $<

