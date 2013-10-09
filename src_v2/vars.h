/* ==========================================================================*/
/*   Copyright (c) 2013       Cullan Howlett & Marc Manera,                  */
/*                            Institute of Cosmology and Gravitation,        */
/*                            University of Portsmouth.                      */
/*                                                                           */
/*   This file is part of PICOLA.                                            */
/*                                                                           */
/*   PICOLA is free software: you can redistribute it and/or modify          */
/*   it under the terms of the GNU General Public License as published by    */
/*   the Free Software Foundation, either version 3 of the License, or       */
/*   (at your option) any later version.                                     */
/*                                                                           */
/*   PICOLA is distributed in the hope that it will be useful,               */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*   GNU General Public License for more details.                            */
/*                                                                           */
/*   You should have received a copy of the GNU General Public License       */
/*   along with PICOLA.  If not, see <http://www.gnu.org/licenses/>.         */
/* ==========================================================================*/

/* ====================================================================================================*/
/* v1.0: This file contains all the external header files required and the global variable definitions.*/
/* v2.0: This file contains extra variables for non-gaussianity and transfer function implementation.  */
/* ====================================================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// GSL libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_integration.h>

// MPI and FFTW libraries
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

// Some definitions
#define  PI          3.14159265358979323846 // PI (obviously)
#define  GRAVITY     6.672e-8               // Newton's constant (in cm^3/g/s^2)
#define  HUBBLE      3.2407789e-18          // Hubble constant (in h/s)

#ifdef SINGLE_PRECISION
typedef float         float_kind;    // Single precision floating types
typedef fftwf_complex complex_kind;  // Single precision complex type
typedef fftwf_plan    plan_kind;     // Single precision FFTW plan
#else
typedef double       float_kind;     // Double precision floating types
typedef fftw_complex complex_kind;   // Double precision complex type
typedef fftw_plan    plan_kind;      // Double precision FFTW plan
#endif

// MPI variables
extern int ierr;             // The return value for mpi routines
extern int NTask;            // The total number of tasks
extern int ThisTask;         // The rank of each task
extern int LeftTask;         // The first neighbouring task on the left containing particles
extern int RightTask;        // The first neighbouring task on the right containing particles
extern MPI_Status status;    // The MPI error status
extern MPI_Request request;  // The continue directive for non-blocking sends

// Global variables for the grids
extern int NTaskWithN;           // The number of tasks that actually have particles
extern int last_slice;           // The last slice of the density/force grids (maybe equal to alloc_local)
extern int * Slab_to_task;       // The task to which each slice is assigned
extern int * Part_to_task;       // The task to which each particle position is assigned
extern int * Local_nx_table;     // The number of slices on each of the tasks
extern int * Local_np_table;     // The number of particle grid slices on each of the tasks
extern float_kind * N11;         // The force grid in the X direction
extern float_kind * N12;         // The force grid in the Y direction
extern float_kind * N13;         // The force grid in the Z direction
extern float_kind * density;     // The density grid
extern ptrdiff_t Local_nx;       // The number of slices on the task
extern ptrdiff_t Local_np;       // The number of particle grid slices on the task
extern ptrdiff_t Total_size;     // The total byte-size of the grids on each processor
extern ptrdiff_t alloc_local;    // The byte-size returned by FFTW required to allocate the density/force grids
extern ptrdiff_t alloc_slice;    // The byte-size of a slice of the density/force grids
extern ptrdiff_t Local_x_start;  // The global start of the slices on the task
extern ptrdiff_t Local_p_start;  // The global start of the particle grid slices on the task
extern complex_kind * P3D;       // Pointer to the complex, FFT'ed density grid (use in-place FFT)
extern complex_kind * FN11;      // Pointer to the complex, FFT'ed N11 force grid (use in-place FFT)
extern complex_kind * FN12;      // Pointer to the complex, FFT'ed N12 force grid (use in-place FFT)
extern complex_kind * FN13;      // Pointer to the complex, FFT'ed N13 force grid (use in-place FFT)
extern plan_kind plan;           // The plan for the in-place FFT of the density grid
extern plan_kind p11,p12,p13;    // Plans for the in-place FFT's of the forces grids 

// Units
extern double G;              // The unit-less Gravitational constant
extern double Hubble;         // The unit-less Hubble constant
extern double UnitMass_in_g;  // The unit mass (in g/cm) used in the code, read in from run parameters
extern double UnitTime_in_s;                   // The unit time (in s) used for the code, calculated from unit length and velocity
extern double UnitLength_in_cm;                // The unit length (in cm/h) used in the code, read in from run parameters file
extern double UnitVelocity_in_cm_per_s;        // The unit velocity (in cm/s) used in the code, read in from run parameters file
extern double InputSpectrum_UnitLength_in_cm;  // The unit length (in cm/h) of the tabulated input spectrum, read in from run parameters

// Cosmological parameters (at z=0)
extern double Fnl;          // The primordial non-gaussianity parameter for local, equilateral or orthogonal
extern double Anorm;        // The normalisation of the power spectrum/ transfer function
extern double Omega;        // The total matter density, CDM+Baryon
extern double Dplus;        // The growth factor
extern double Sigma8;       // The normalisation of the power spectrum 
extern double FnlTime;      // The scale factor at which fnl kicks in
extern double InitTime;     // The scale factor at the final redshift
extern double DstartFnl;    // The growth factor for the initial potential 
extern double ShapeGamma;   // The paramerisation of the Efstathiou power spectrum
extern double OmegaBaryon;  // The baryonic matter density
extern double HubbleParam;  // The normalised Hubble parameter, h=H/100
extern double OmegaLambda;         // The Cosmological Constant
extern double Fnl_Redshift;        // The redshift at which the nongaussian f_nl potential is computed
extern double Init_Redshift;       // The redshift at which to begin timestepping
extern double Final_Redshift;      // The final redshift of the simulation
extern double PrimordialIndex;     // The spectral index, n_s

#ifdef GENERIC_FNL
// Kernels to include general non-gaussian models
extern int NKernelTable;           // The length of the kernel lookup table
extern struct kern_table {         // The kernel lookup table
  double Coef;
  double ker0;
  double kerA;
  double kerB;
} *KernelTable;
#endif

// Particle data and pointers
#ifdef MEMORY_MODE
extern float * Disp[3];    // Vectors to hold the particle displacements each timestep
extern float * ZA[3];      // Vectors to hold the Zeldovich displacements before particle initialisation
extern float * LPT[3];     // Vectors to hold the 2LPT displacements before particle initialisation
extern struct part_data { 
  unsigned long long ID;   // The Particle ID
  float Dz[3];             // The Zeldovich displacment of the particle in the X, Y and Z directions
  float D2[3];             // The 2LPT displacment of the particle in the X, Y and Z directions
  float Pos[3];            // The position of the particle in the X, Y and Z directions
  float Vel[3];            // The velocity of the particle in the X, Y and Z directions
} *P;
#else
extern float_kind * Disp[3];  // vectors to hold the particle displacements each timestep
extern float_kind * ZA[3];    // Vectors to hold the Zeldovich displacements before particle initialisation
extern float_kind * LPT[3];   // Vectors to hold the 2LPT displacements before particle initialisation
extern struct part_data {
  unsigned long long ID;      // The particle ID
  float_kind Dz[3];           // The Zeldovich displacment of the particle in the X, Y and Z directions
  float_kind D2[3];           // The 2LPT displacment of the particle in the X, Y and Z directions
  float_kind Pos[3];          // The position of the particle in the X, Y and Z directions
  float_kind Vel[3];          // The velocity of the particle in the X, Y and Z directions
} *P;
#endif

extern struct io_header_1
{
  unsigned int npart[6];       // npart[1] gives the number of particles in the present file, other particle types are ignored
  double mass[6];              // mass[1] gives the particle mass
  double time;                 // time (=cosmological scale factor) of snapshot 
  double redshift;             // redshift of snapshot
  int flag_sfr;                // flags whether star formation is used (not available in L-Gadget2)
  int flag_feedback;           // flags whether feedback from star formation is included
  unsigned int npartTotal[6];  // npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                               // the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder.
  int flag_cooling;            // flags whether radiative cooling is included
  int num_files;               // determines the number of files that are used for a snapshot
  double BoxSize;              // Simulation box size (in code units)
  double Omega0;               // matter density
  double OmegaLambda;          // vacuum energy density
  double HubbleParam;          // little 'h'
  int flag_stellarage;         // flags whether the age of newly formed stars is recorded and saved
  int flag_metals;             // flags whether metal enrichment is included
  int hashtabsize;             // gives the size of the hashtable belonging to this snapshot file
  char fill[84];	       // fills to 256 Bytes
} header, header1;

// Simulation variables
// ====================
extern char FileBase[100];   // The base output filename
extern char OutputDir[100];  // The output directory
extern int Nmesh;            // The size of the displacement, density and force grids (in 1-D) 
extern int Nsample;          // The number of particles (in 1-D)
extern int UseCOLA;          // Whether or not to use the COLA modifications
extern int NumFilesWrittenInParallel;  // The maximum number of files to be written out in parallel
extern unsigned int NumPart;           // The number of particles on each processor
extern unsigned long long TotNumPart;  // The total number of particles in the simulation
extern double Box;                     // The edge length of the simulation
extern double Buffer;                  // The amount of extra memory of each processor to compensate for moving particles

// 2LPT specific
extern char FileWithInputSpectrum[500];  // The file containing the input power spectrum
extern char FileWithInputTransfer[500];  // The file containing the input transfer function
extern char FileWithInputKernel[500];    // The file containing the input nongaussian kernel
extern int Seed;                         // The random seed to generate to realisation
extern int SphereMode;       // Whether to use a sphere or a cube in k-space
extern int WhichSpectrum;    // Which power spectrum to use
extern int WhichTransfer;    // Which transfer function to use

// COLA specific
extern int fullT;    // The time dependence of the velocity, hardcoded (see README)
extern int StdDA;    // The time dependence of the displacement, hardcoded (see README)
extern int filter;   // Whether or not to apply a smoothing filter to the forces, hardcoded
extern int nsteps;   // The number of timesteps to perform
extern double nLPT;   // Parameterisation of the time dependence of the velocity, hardcoded (see README)
extern double Scale;  // The smoothing scale for the smoothing filter, hardcoded
