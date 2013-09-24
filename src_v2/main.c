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

/* ======================================================*/
/* This file contains the main driver routine for PICOLA.*/
/* ======================================================*/
     
#include "vars.h"
#include "proto.h"

int main(int argc, char **argv) {
   
  // Set up MPI
  // ==========
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#ifdef SINGLE_PRECISION
  fftwf_mpi_init();
#else
  fftw_mpi_init();
#endif

  if(argc < 2) {
    if(ThisTask == 0) {
      fprintf(stdout, "Input parameters not found\n");
      fprintf(stdout, "Call with <ParameterFile>\n");
    }
    ierr = MPI_Finalize();
    exit(0);
  }
   
  // Read the run parameters and setup code
  // ======================================
  int stepDistr;   
  int subtractLPT;
  double da=0;

  read_parameterfile(argv[1]);
  
  if (UseCOLA == 1){
    subtractLPT = 1; 
    stepDistr   = 0;
    StdDA       = 0;
  } else{
    subtractLPT = 0; 
    stepDistr   = 1;
    StdDA       = 2;
  }
  if (StdDA == 0){
    fullT = 1;
    nLPT  = -2.5;
  }
  filter = 0;              // Whether or not to smooth the forces
  Scale  = 2.*M_PI/Box;    // The force smoothing scale 

  if(ThisTask == 0) {
    printf("Run Parameters\n");
    printf("==============\n");
    printf("Cosmology:\n");
    printf("  Omega Matter(z=0) = %lf\n",Omega);
    printf("  Omega Baryon(z=0) = %lf\n",OmegaBaryon);
    printf("  Hubble Parameter(z=0) = %lf\n",HubbleParam);
    printf("  Sigma8(z=0) = %lf\n",Sigma8);
#ifndef GAUSSIAN
    printf("  F_nl = %lf\n",Fnl);
#endif
    printf("  Primordial Index = %lf\n",PrimordialIndex);
    printf("  Initial Redshift  = %lf\n",Init_Redshift);
    printf("  Final Redshift    = %lf\n",Final_Redshift);
#ifndef GAUSSIAN
    printf("  F_nl Redshift  = %lf\n",Fnl_Redshift);
#endif
    printf("Simulation:\n");
    printf("  Nmesh = %d\n", Nmesh);
    printf("  Nsample = %d\n", Nsample);
    printf("  Boxsize = %lf\n", Box);
    printf("  Buffer Size = %lf\n", Buffer);
    switch(WhichSpectrum) {
      case 0:
        switch (WhichTransfer) {
          case 1:
            printf("  Using Eisenstein & Hu Transfer Function\n");
            break;
          case 2:
            printf("  Using Tabulated Transfer Function\n");
            break;
          default:
            printf("  Using Efstathiou Transfer Function\n");
            break;
        }
        break;
      case 1:
        printf("  Using Eisenstein & Hu Power Spectrum\n");
        break;
      case 2:
        printf("  Using Tabulated Power Spectrum\n");
        break;   
      default:
        printf("  Using Efstathiou Power Spectrum\n");
        break;
    }      
    printf("  Number of Timesteps = %d\n",nsteps);
    if (UseCOLA) {
      printf("  Using COLA method\n\n");
    } else {
      printf("  Using Standard PM method\n\n");
    }
    fflush(stdout);
  }   
  
  // Initial and final scale factors:
  double ai=1.0/(1.0+Init_Redshift);
  double af=1.0/(1.0+Final_Redshift);
    
  if (stepDistr == 0) da=(af-ai)/((double)nsteps);
  if (stepDistr == 1) da=(log(af)-log(ai))/((double)nsteps);
  if (stepDistr == 2) da=(CosmoTime(af)-CosmoTime(ai))/((double)nsteps);

  set_units();

  if (ThisTask == 0) {
    printf("Initialising Transfer Function/Power Spectrum\n");
    printf("=============================================\n");
  }
  initialize_transferfunction();
  initialize_powerspectrum();
  initialize_ffts();
  initialize_parts();

  if(ThisTask == 0) {
    printf("Creating initial conditions\n");
    printf("===========================\n");
    fflush(stdout);
  }

  // Create the calculate the Zeldovich and 2LPT displacements and create the initial conditions
  // ===========================================================================================
  int i, j, k, m;
  unsigned int n, coord;
  double A=ai;                // This is the scale factor which we'll be advancing below.
  double Di=growthD(1.0, A);  // initial growth factor
  double Di2=growthD2(A);     // initial 2nd order growth factor  
  double Dv=DprimeQ(A,1.0);   // T[D_{za}]=dD_{za}/dy
  double Dv2=growthD2v(A);    // T[D_{2lpt}]=dD_{2lpt}/dy

  displacement_fields();
    
  P = (struct part_data *) malloc((int)(ceil(NumPart*Buffer))*sizeof(struct part_data));

  // Generate the initial particle positions and velocities
  // If subtractLPT = 0 (non-COLA), then velocity is ds/dy, which is simply the 2LPT IC.
  // Else set vel = 0 if we subtract LPT. This is the same as the action of the operator L_- from TZE, as initial velocities are in 2LPT.
  for(i=0; i<Local_np; i++) {
    for (j=0; j<Nsample; j++) {
      for (k=0; k<Nsample; k++) {
        coord = (i * Nsample + j) * Nsample + k;
           
        for (m=0; m<3; m++) {
          P[coord].Dz[m] = ZA[m][coord];
          P[coord].D2[m] = LPT[m][coord];
          if (subtractLPT == 0) {
            P[coord].Vel[m]=P[coord].Dz[m]*Dv+P[coord].D2[m]*Dv2;
          } else {
            P[coord].Vel[m] = 0.0;
          }
        }

        P[coord].Pos[0] = periodic_wrap((i+Local_p_start)*(Box/Nsample)+P[coord].Dz[0]*Di+P[coord].D2[0]*Di2);
        P[coord].Pos[1] = periodic_wrap(j*(Box/Nsample)+P[coord].Dz[1]*Di+P[coord].D2[1]*Di2);
        P[coord].Pos[2] = periodic_wrap(k*(Box/Nsample)+P[coord].Dz[2]*Di+P[coord].D2[2]*Di2);
      }
    }
  }

  for (i=0; i<3; i++) {
    free(ZA[i]);
    free(LPT[i]);
  }

  // Now, we get to the N-Body part where we evolve with time via the Kick-Drift-Kick Method
  // =======================================================================================
  int timeStep;
  double AF=0,AI,AC,AFF=0;
  double growth1   = Di;
  double growth1L2 = Di2;

  // The density grid and force grids  and associated fftw plans
#ifndef MEMORY_MODE
  density = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N11  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N12  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N13  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  P3D  = (complex_kind*)density;
  FN11 = (complex_kind*)N11;
  FN12 = (complex_kind*)N12;
  FN13 = (complex_kind*)N13;
#ifdef SINGLE_PRECISION
  plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p11  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p12  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p13  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else
  plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p11  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p12  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p13  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#endif
#endif

 if(ThisTask == 0) {
    printf("Beginning timestepping\n");
    printf("======================\n");
    fflush(stdout);
  }
  
  // AI stores the scale factor to which the velocities have been kicked to. Initially it's just A.
  AI=A;
  for (timeStep=0;timeStep<=nsteps;timeStep++){
    
    // AFF is the scale factor to which we should drift the particle positions.
    // AF is the scale factor to which we should kick the particle velocities.
    if (stepDistr == 0) AFF=A+da;
    if (stepDistr == 1) AFF=A*exp(da);
    if (stepDistr == 2) AFF=AofTime(CosmoTime(A)+da);

    // half time-step for final kick
    if (timeStep == nsteps) {
      AF=A; 
    } else { 
      // Set to mid-point of interval. In the infinitesimal timestep limit, these choices are identical. 
      // How one chooses the mid-point when not in that limit is really an extra degree of freedom in the code 
      // but Tassev et al. report negligible effects from the different choices below. 
      // Hence, this is not exported as an extra switch at this point.
      if (stepDistr == 0) AF=A+da*0.5;
      if (stepDistr == 1) AF=A*exp(da*0.5);
      if (stepDistr == 2) AF=AofTime((CosmoTime(AFF)+CosmoTime(A))*0.5); 
    }
    
    if (ThisTask == 0) {
      printf("Iteration = %d\n------------------\n",timeStep+1);
      printf("a = %lf\n",A);
      printf("z = %lf\n",1.0/A-1.0);
      fflush(stdout);
    }

    // First we check whether all the particles are on the correct processor after the last time step/
    // original 2LPT displacement and move them if not
    if (ThisTask == 0) printf("Moving particles across task boundaries...\n");
    MoveParticles();

#ifdef MEMORY_MODE
    density = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
    P3D  = (complex_kind*)density;
#ifdef SINGLE_PRECISION
    plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else
    plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
#endif
#endif

    // Then we do the Cloud-in-Cell assignment to get the density grid and FFT it.  
    if (ThisTask == 0) printf("Calculating density using Cloud-in-Cell...\n");
    PtoMesh();

#ifdef MEMORY_MODE
    N11  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
    N12  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
    N13  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
    FN11 = (complex_kind*)N11;
    FN12 = (complex_kind*)N12;
    FN13 = (complex_kind*)N13;
#ifdef SINGLE_PRECISION
    p11  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
    p12  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
    p13  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else
    p11  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
    p12  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
    p13  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#endif
#endif
    
    // This returns N11,N12,N13 which hold the components of
    // the vector (grad grad^{-2} density) on a grid.
    if (ThisTask == 0) printf("Calculating forces...\n");
    Forces(); 

#ifdef MEMORY_MODE
    free(density);
    for (i=0; i<3; i++) Disp[i] = (float *)malloc(NumPart*sizeof(float));
#ifdef SINGLE_PRECISION
    fftwf_destroy_plan(plan);
#else 
    fftw_destroy_plan(plan);
#endif
#else
    for (i=0; i<3; i++) Disp[i] = (float_kind *)malloc(NumPart*sizeof(float_kind));
#endif
    
    // Now find the accelerations at the particle positions using 3-linear interpolation. 
    if (ThisTask == 0) printf("Calculating accelerations...\n");
    MtoParticles();

#ifdef MEMORY_MODE
  free(N11);
  free(N12);
  free(N13);  
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(p11);
  fftwf_destroy_plan(p12);
  fftwf_destroy_plan(p13);
#else
  fftw_destroy_plan(p11);
  fftw_destroy_plan(p12);
  fftw_destroy_plan(p13);
#endif
#endif
    
    // Calculate the mean displacement and subtract later.
    if (ThisTask == 0) printf("Calculating mean of displacements...\n");
    double sumDx=0,sumDy=0,sumDz=0;
    for(n=0; n<NumPart; n++) {
      sumDx += Disp[0][n];
      sumDy += Disp[1][n];
      sumDz += Disp[2][n];
    }

    // Make sumDx, sumDy and sumDz global averages
    ierr = MPI_Allreduce(MPI_IN_PLACE,&sumDx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(MPI_IN_PLACE,&sumDy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(MPI_IN_PLACE,&sumDz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  
    
    sumDx /= (double)TotNumPart; // We will subtract these below to conserve momentum. 
    sumDy /= (double)TotNumPart;
    sumDz /= (double)TotNumPart; 

    if (ThisTask == 0) {
      printf("Kicking the particles...\n");
      fflush(stdout);
    }

    // Kick
    // ===============
    double dda;
    double q1,q2;
    double ax,ay,az;
    double sumx=0,sumy=0,sumz=0; 
    double Om143=pow(Omega/(Omega+(1-Omega)*A*A*A),1./143.);
    
    if (StdDA == 0) {
      dda=Sphi(AI,AF,A);
    } else if (StdDA == 1) {
      dda=(AF-AI)*A/Qfactor(A);
    } else {
      dda=SphiStd(AI,AF);
    }  
    
    q2=1.5*Omega*growth1*growth1*(1.0+7./3.*Om143)*A; // T^2[D_{2lpt}]=d^2 D_{2lpt}/dy^2
    q1=1.5*Omega*growth1*A;                           // T^2[D_{ZA}]=d^2 D_{ZA}/dy^2
    
    for(n=0; n<NumPart; n++) {

      Disp[0][n] -= sumDx;
      Disp[1][n] -= sumDy;
      Disp[2][n] -= sumDz;

      ax=-1.5*Omega*Disp[0][n]-subtractLPT*(P[n].Dz[0]*q1+P[n].D2[0]*q2)/A;
      ay=-1.5*Omega*Disp[1][n]-subtractLPT*(P[n].Dz[1]*q1+P[n].D2[1]*q2)/A;
      az=-1.5*Omega*Disp[2][n]-subtractLPT*(P[n].Dz[2]*q1+P[n].D2[2]*q2)/A;

      P[n].Vel[0] += ax*dda;
      P[n].Vel[1] += ay*dda;
      P[n].Vel[2] += az*dda;

      sumx += P[n].Vel[0];
      sumy += P[n].Vel[1];
      sumz += P[n].Vel[2];
    }

    for (i=0; i<3; i++) free(Disp[i]);

    // Make sumx, sumy and sumz global averages
    ierr = MPI_Allreduce(MPI_IN_PLACE,&sumx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(MPI_IN_PLACE,&sumy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    ierr = MPI_Allreduce(MPI_IN_PLACE,&sumz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  
    
    sumx /= (double)TotNumPart;  // We will subtract these below to conserve momentum. 
    sumy /= (double)TotNumPart;  // Should be conserved, but just in case 3-linear interpolation makes a problem.
    sumz /= (double)TotNumPart;  // Never checked whether this makes a difference.

    if (timeStep == nsteps) {

      if (ThisTask == 0) {
        printf("Iteration %d finished\n------------------\n\n", timeStep+1);
        printf("Timestepping finished\n\n");
        fflush(stdout);
      }
    
      // At final timestep, add back LPT velocities if we had subtracted them. 
      // This corresponds to L_+ operator in TZE.
      Dv  = DprimeQ(A,1.0);  // dD_{za}/dy
      Dv2 = growthD2v(A);    // dD_{2lpt}/dy

      for(n=0; n<NumPart; n++) {
        P[n].Vel[0] += -sumx+(P[n].Dz[0]*Dv+P[n].D2[0]*Dv2)*subtractLPT;
        P[n].Vel[1] += -sumy+(P[n].Dz[1]*Dv+P[n].D2[1]*Dv2)*subtractLPT;
        P[n].Vel[2] += -sumz+(P[n].Dz[2]*Dv+P[n].D2[2]*Dv2)*subtractLPT;
      }

      goto finalize; // Sorry for "goto" :)
    }
    
    if (ThisTask == 0) {
      printf("Drifting the particles...\n");
      fflush(stdout);
    }

    // Drift
    // =============
    double dyyy;
    double da1,da2;

    AC = AF;
    AF = AFF;
    
    if (StdDA == 0) {
      dyyy=Sq(A,AF,AC);
    } else if (StdDA == 1) {
      dyyy=(AF-A)/Qfactor(AC);
    } else {
      dyyy=SqStd(A,AF);
    }

    da1=growthD(1.0, AF)-growth1;    // change in D
    da2=growthD2(AF)-growth1L2; // change in D_{2lpt}
    
    for(n=0; n<NumPart; n++) {
        P[n].Pos[0] += (P[n].Vel[0]-sumx)*dyyy;
        P[n].Pos[1] += (P[n].Vel[1]-sumy)*dyyy;
        P[n].Pos[2] += (P[n].Vel[2]-sumz)*dyyy;

        P[n].Pos[0] = periodic_wrap(P[n].Pos[0]+subtractLPT*(P[n].Dz[0]*da1+P[n].D2[0]*da2));
        P[n].Pos[1] = periodic_wrap(P[n].Pos[1]+subtractLPT*(P[n].Dz[1]*da1+P[n].D2[1]*da2));
        P[n].Pos[2] = periodic_wrap(P[n].Pos[2]+subtractLPT*(P[n].Dz[2]*da1+P[n].D2[2]*da2));
    }

    // Step in time
    // ================
    A  = AF;   // WRT to the above name change, A  = AFF
    AI = AC;   // WRT to the above name change, AI = AF

    growth1   = growthD(1.0, A);
    growth1L2 = growthD2(A);

    if (ThisTask == 0) {
      printf("Iteration %d finished\n------------------\n\n", timeStep+1);
      fflush(stdout);
    }
     
    ierr = MPI_Barrier(MPI_COMM_WORLD);

  }

  // Here is the last little bit
  // ===========================
  finalize:

  if (ThisTask == 0) {
    printf("Finishing up\n");
    printf("============\n");
    fflush(stdout);
  }
    
  // Now convert velocities to v_{rsd}\equiv (ds/d\eta)/(a H(a))
  velRSD(A);
    
  // Output a slice just for the sake of doing something with P.
  if (ThisTask == 0) {
    printf("Converting to RSD velocities...\n");
    printf("Outputting particles...\n"); 
  }
  slice();
  print_spec();
  fflush(stdout);

  free_powertable();
  free_transfertable();
#ifdef GENERIC_FNL
  free(KernelTable);
#endif

  free(P);
  free(Slab_to_task);
  free(Part_to_task);
  free(Local_nx_table);
  free(Local_np_table);
#ifndef MEMORY_MODE
  free(density);
  free(N11);
  free(N12);
  free(N13);  
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(plan);
  fftwf_destroy_plan(p11);
  fftwf_destroy_plan(p12);
  fftwf_destroy_plan(p13);
#else
  fftw_destroy_plan(plan);
  fftw_destroy_plan(p11);
  fftw_destroy_plan(p12);
  fftw_destroy_plan(p13);
#endif
#endif

#ifdef SINGLE_PRECISION
  fftwf_mpi_cleanup();
#else
  fftw_mpi_cleanup();
#endif

  if (ThisTask == 0) printf("Done :)\n");

  MPI_Finalize();   

  return 0;
}
