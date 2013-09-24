CC        = mpiCC
OPTIMIZE  = -O3 -Wall

# Precision
#==========

#SINGLE_PRECISION = -DSINGLE_PRECISION	# Single precision floats and FFTW (else use double precision)
#OPTIONS += $(SINGLE_PRECISION)

MEMORY_MODE = -DMEMORY_MODE		# Save memory by making sure to allocate and deallocate arrays only when we need them
OPTIONS += $(MEMORY_MODE)		# and by making the particle data single precision

# Gaussian/Non-gaussian initial conditions
# ========================================

GAUSSIAN = -DGAUSSIAN          # Switch this if you want gaussian initial conditions (fnl otherwise)
OPTIONS += $(GAUSSIAN) 

#LOCAL_FNL = -DLOCAL_FNL        # Switch this if you want only local non-gaussianities
#OPTIONS += $(LOCAL_FNL)        # NOTE this option is only for invariant inital power spectrum
                                # for local with ns != 1 use DGENERIC_FNL and input_kernel_local.txt

#EQUIL_FNL = -DEQUIL_FNL        # Switch this if you want equilateral Fnl
#OPTIONS += $(EQUIL_FNL)        # NOTE this option is only for invariant inital power spectrum
                                # for local with ns != 1 use DGENERIC_FNL and input_kernel_equil.txt

#ORTHO_FNL = -DORTHO_FNL        # Switch this if you want ortogonal Fnl
#OPTIONS += $(ORTHO_FNL)        # NOTE this option is only for invariant inital power spectrum
                                # for local with ns != 1 use DGENERIC_FNL and input_kernel_ortog.txt

#GENERIC_FNL += -DGENERIC_FNL   # Switch this if you want generic Fnl implementation
#OPTIONS += $(GENERIC_FNL)      # This option allows for ns != 1 and should include an input_kernel_file.txt 
                                # containing the coefficients for the generic kernel 
                                # see README and Manera et al astroph/NNNN.NNNN
                                # For local, equilateral and orthogonal models you can use the provided files
                                # input_kernel_local.txt, input_kernel_equil.txt, input_kernel_orthog.txt 

# Run some checks on option compatability
# =======================================
ifdef GAUSSIAN
ifdef LOCAL_FNL
  $(error GAUSSIAN AND LOCAL_FNL are not compatible, change Makefile)
endif
ifdef EQUIL_FNL
  $(error GAUSSIAN AND EQUIL_FNL are not compatible, change Makefile)
endif
ifdef ORTHO_FNL
  $(error GAUSSIAN AND ORTHO_FNL are not compatible, change Makefile)
endif
else
ifndef LOCAL_FNL 
ifndef EQUIL_FNL
ifndef ORTHO_FNL 
ifndef GENERIC_FNL
  $(error if not using GAUSSIAN then must select some type of non-gaussianity (LOCAL_FNL, EQUIL_FNL, ORTHO_FNL, GENERIC_FNL), change Makefile)
endif
endif
endif
endif
endif

ifdef GENERIC_FNL 
ifdef LOCAL_FNL 
   $(error GENERIC_FNL AND LOCAL_FNL are not compatible, choose one in Makefile) 
endif 
ifdef EQUIL_FNL 
   $(error GENERIC_FNL AND EQUIL_FNL are not compatible, choose one in Makefile) 
endif 
ifdef ORTHO_FNL 
   $(error GENERIC_FNL AND ORTHO_FNL are not compatible, choose one in Makefile) 
endif 
endif 

ifdef LOCAL_FNL
ifdef EQUIL_FNL
   $(error LOCAL_FNL AND EQUIL_FNL are not compatible, choose one or the other in Makefile) 
endif
ifdef ORTHO_FNL
   $(error LOCAL_FNL AND ORTHO_FNL are not compatible, choose one or the other in Makefile) 
endif
endif

ifdef EQUIL_FNL
ifdef ORTHO_FNL
   $(error EQUIL_FNL AND ORTHO_FNL are not compatible, choose one or the other in Makefile) 
endif
endif

# Setup libraries and compile the code
# ====================================
ifdef SINGLE_PRECISION
FFTW_INCL = -I/opt/gridware/libs/gcc/fftw3/3_3_2/include/
FFTW_LIBS = -L/opt/gridware/libs/gcc/fftw3/3_3_2/lib/ -lfftw3f_mpi -lfftw3f
else
FFTW_INCL = -I/opt/gridware/libs/gcc/fftw3/3.3.3/include/
FFTW_LIBS = -L/opt/gridware/libs/gcc/fftw3/3.3.3/lib/ -lfftw3_mpi -lfftw3
endif

GSL_INCL  = -I/opt/gridware/libs/gcc/gsl/1.14/include/
GSL_LIBS  = -L/opt/gridware/libs/gcc/gsl/1.14/lib/  -lgsl -lgslcblas

MPILIBS = -L/opt/gridware/mpi/gcc/openmpi/1_4_3/lib/ -lmpi

LIBS   =   -lm $(MPILIBs) $(FFTW_LIBS) $(GSL_LIBS)

CFLAGS =   $(OPTIMIZE) $(FFTW_INCL) $(GSL_INCL) $(OPTIONS)

OBJS   = src_v2/main.o src_v2/cosmo.o src_v2/auxPM.o src_v2/2LPT.o src_v2/power.o src_v2/vars.o src_v2/read_param.o
ifdef GENERIC_FNL
OBJS += src_v2/kernel.o
endif

INCL   = src_v2/vars.h src_v2/proto.h  Makefile

all: $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o  PICOLA

$(OBJS): $(INCL) 

clean:
	rm -f src_v2/*.o src_v2/*~ files/*~ *~ PICOLA 
