MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.mld2p4
MLDMODDIR=$(MLDDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDMODDIR) $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

LINKOPT=
EXEDIR=./runs

all: mld_s_pde3d mld_d_pde3d mld_s_pde2d mld_d_pde2d 

mld_d_pde3d: mld_d_pde3d.o  data_input.o
	$(FLINK) $(LINKOPT) mld_d_pde3d.o  data_input.o -o mld_d_pde3d $(MLD_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_d_pde3d $(EXEDIR)

mld_s_pde3d: mld_s_pde3d.o   data_input.o
	$(FLINK) $(LINKOPT) mld_s_pde3d.o  data_input.o -o mld_s_pde3d $(MLD_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_s_pde3d $(EXEDIR)

mld_d_pde2d: mld_d_pde2d.o  data_input.o
	$(FLINK) $(LINKOPT) mld_d_pde2d.o  data_input.o -o mld_d_pde2d $(MLD_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_d_pde2d $(EXEDIR)

mld_s_pde2d: mld_s_pde2d.o   data_input.o
	$(FLINK) $(LINKOPT) mld_s_pde2d.o  data_input.o -o mld_s_pde2d $(MLD_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_s_pde2d $(EXEDIR)

mld_d_pde3d.o mld_s_pde3d.o mld_d_pde2d.o mld_s_pde2d.o: data_input.o 

check: all
	cd runs && ./mld_d_pde2d <mld_pde2d.inp && ./mld_s_pde2d<mld_pde2d.inp


clean: 
	/bin/rm -f data_input.o mld_d_pde3d.o mld_s_pde3d.o mld_d_pde2d.o mld_s_pde2d.o *$(.mod)\
        $(EXEDIR)/mld_d_pde3d  $(EXEDIR)/mld_s_pde3d $(EXEDIR)/mld_d_pde2d  $(EXEDIR)/mld_s_pde2d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



