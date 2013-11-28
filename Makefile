# The executable name
# ===================
EXEC    = PICOLA

# Choose the machine you are running on. Currently versions for SCIAMA and DARWIN are implemented
# ===============================================================================================
MACHINE = SCIAMA
#MACHINE = DARWIN

# Options for optimization
# ========================
OPTIMIZE  = -O3 -Wall

# Various C preprocessor directives that change the way PICOLA is made
# ====================================================================

#SINGLE_PRECISION = -DSINGLE_PRECISION	# Single precision floats and FFTW (else use double precision)
#OPTIONS += $(SINGLE_PRECISION)

MEMORY_MODE = -DMEMORY_MODE		# Save memory by making sure to allocate and deallocate arrays only when we need them
OPTIONS += $(MEMORY_MODE)		# and by making the particle data single precision

#LIGHTCONE = -DLIGHTCONE                 # Builds a lightcone based on the run parameters and only outputs particles
#OPTIONS += $(LIGHTCONE)                 # at a given timestep if they have entered the lightcone 

GAUSSIAN = -DGAUSSIAN                   # Switch this if you want gaussian initial conditions (fnl otherwise)
OPTIONS += $(GAUSSIAN) 

#LOCAL_FNL = -DLOCAL_FNL                # Switch this if you want only local non-gaussianities
#OPTIONS += $(LOCAL_FNL)                # NOTE this option is only for invariant inital power spectrum
                                        # for local with ns != 1 use DGENERIC_FNL and input_kernel_local.txt

#EQUIL_FNL = -DEQUIL_FNL                # Switch this if you want equilateral Fnl
#OPTIONS += $(EQUIL_FNL)                # NOTE this option is only for invariant inital power spectrum
                                        # for local with ns != 1 use DGENERIC_FNL and input_kernel_equil.txt

#ORTHO_FNL = -DORTHO_FNL                # Switch this if you want ortogonal Fnl
#OPTIONS += $(ORTHO_FNL)                # NOTE this option is only for invariant inital power spectrum
                                        # for local with ns != 1 use DGENERIC_FNL and input_kernel_ortog.txt

#GENERIC_FNL += -DGENERIC_FNL           # Switch this if you want generic Fnl implementation
#OPTIONS += $(GENERIC_FNL)              # This option allows for ns != 1 and should include an input_kernel_file.txt 
                                        # containing the coefficients for the generic kernel 
                                        # see README and Manera et al astroph/NNNN.NNNN
                                        # For local, equilateral and orthogonal models you can use the provided files
                                        # input_kernel_local.txt, input_kernel_equil.txt, input_kernel_orthog.txt 

GADGET_STYLE = -DGADGET_STYLE          # Writes all the output in Gadget's '1' style format, with the corresponding
OPTIONS += $(GADGET_STYLE)             # header and correct velocity units


# Nothing below here should need changing unless you are adding in/modifying libraries for existing or new machines
# =================================================================================================================

# Run some checks on option compatability
# =======================================
ifdef GAUSSIAN
ifdef LOCAL_FNL
  $(error ERROR: GAUSSIAN AND LOCAL_FNL are not compatible, change Makefile)
endif
ifdef EQUIL_FNL
  $(error ERROR: GAUSSIAN AND EQUIL_FNL are not compatible, change Makefile)
endif
ifdef ORTHO_FNL
  $(error ERROR: GAUSSIAN AND ORTHO_FNL are not compatible, change Makefile)
endif
else
ifndef LOCAL_FNL 
ifndef EQUIL_FNL
ifndef ORTHO_FNL 
ifndef GENERIC_FNL
  $(error ERROR: if not using GAUSSIAN then must select some type of non-gaussianity (LOCAL_FNL, EQUIL_FNL, ORTHO_FNL, GENERIC_FNL), change Makefile)
endif
endif
endif
endif
endif

ifdef GENERIC_FNL 
ifdef LOCAL_FNL 
   $(error ERROR: GENERIC_FNL AND LOCAL_FNL are not compatible, choose one in Makefile) 
endif 
ifdef EQUIL_FNL 
   $(error ERROR: GENERIC_FNL AND EQUIL_FNL are not compatible, choose one in Makefile) 
endif 
ifdef ORTHO_FNL 
   $(error ERROR: GENERIC_FNL AND ORTHO_FNL are not compatible, choose one in Makefile) 
endif 
endif 

ifdef LOCAL_FNL
ifdef EQUIL_FNL
   $(error ERROR: LOCAL_FNL AND EQUIL_FNL are not compatible, choose one or the other in Makefile) 
endif
ifdef ORTHO_FNL
   $(error ERROR: LOCAL_FNL AND ORTHO_FNL are not compatible, choose one or the other in Makefile) 
endif
endif

ifdef EQUIL_FNL
ifdef ORTHO_FNL
   $(error ERROR: EQUIL_FNL AND ORTHO_FNL are not compatible, choose one or the other in Makefile) 
endif
endif

ifdef LIGHTCONE
ifdef GADGET_STYLE
   $(warning WARNING: GADGET_STYLE outputting is not supported for lightcone simulations, will revert to standard ASCII format)
endif
endif

# Setup libraries and compile the code
# ====================================
ifeq ($(MACHINE),SCIAMA)
  CC = mpiCC
ifdef SINGLE_PRECISION
  FFTW_INCL = -I/opt/gridware/libs/gcc/fftw3/3_3_2/include/
  FFTW_LIBS = -L/opt/gridware/libs/gcc/fftw3/3_3_2/lib/ -lfftw3f_mpi -lfftw3f
else
  FFTW_INCL = -I/opt/gridware/libs/gcc/fftw3/3.3.3/include/
  FFTW_LIBS = -L/opt/gridware/libs/gcc/fftw3/3.3.3/lib/ -lfftw3_mpi -lfftw3
endif
  GSL_INCL  = -I/opt/gridware/libs/gcc/gsl/1.14/include/
  GSL_LIBS  = -L/opt/gridware/libs/gcc/gsl/1.14/lib/  -lgsl -lgslcblas
  MPI_INCL  = -I/opt/gridware/mpi/gcc/openmpi/1_4_3/include
  MPI_LIBS  = -L/opt/gridware/mpi/gcc/openmpi/1_4_3/lib/ -lmpi
endif

ifeq ($(MACHINE),DARWIN)
  CC = mpiicc	
ifdef SINGLE_PRECISION
  FFTW_INCL = -I/usr/local/Cluster-Apps.sandybridge/fftw/intel/3.3.3/include
  FFTW_LIBS = -L/usr/local/Cluster-Apps.sandybridge/fftw/intel/3.3.3/lib -lfftw3f_mpi -lfftw3f
else
  FFTW_INCL = -I/usr/local/Cluster-Apps.sandybridge/fftw/intel/3.3.3/include
  FFTW_LIBS = -L/usr/local/Cluster-Apps.sandybridge/fftw/intel/3.3.3/lib -lfftw3_mpi -lfftw3
endif
  GSL_INCL  = -I/usr/local/Cluster-Apps/gsl/1.9/include/
  GSL_LIBS  = -L/usr/local/Cluster-Apps/gsl/1.9/lib/  -lgsl -lgslcblas
  MPI_INCL  = -L/usr/local/Cluster-Apps/intel/impi/3.1/include
  MPI_LIBS  = -L/usr/local/Cluster-Apps/intel/impi/3.1/lib -lmpi
endif

LIBS   =   -lm $(MPI_LIBs) $(FFTW_LIBS) $(GSL_LIBS)

CFLAGS =   $(OPTIMIZE) $(FFTW_INCL) $(GSL_INCL) $(MPI_INCL) $(OPTIONS)

OBJS   = src_v3/main.o src_v3/cosmo.o src_v3/auxPM.o src_v3/2LPT.o src_v3/power.o src_v3/vars.o src_v3/read_param.o
ifdef GENERIC_FNL
OBJS += src_v3/kernel.o
endif
ifdef LIGHTCONE
OBJS += src_v3/lightcone.o
endif

INCL   = src_v3/vars.h src_v3/proto.h  Makefile

all: $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f src_v3/*.o src_v3/*~ *~ $(EXEC)
