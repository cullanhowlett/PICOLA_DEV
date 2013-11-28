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
    exit(1);
  }
   
  // Read the run parameters and setup code
  // ======================================
  int i, stepDistr;

  read_parameterfile(argv[1]);
  read_outputs();
  set_units();
#ifdef LIGHTCONE
  set_lightcone();
#endif
  
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
#ifndef GAUSSIAN
    printf("  F_nl Redshift  = %lf\n",Fnl_Redshift);
#endif
    printf("\nSimulation:\n");
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
    if (UseCOLA) {
      printf("  Using COLA method\n");
    } else {
      printf("  Using Standard PM method\n");
    }
#ifdef LIGHTCONE
    printf("\nLightcone:\n");
    printf("  Maximum Comoving Radius = %lf\n", Light/Hubble*SphiStd(1.0/(1.0+OutputList[0].Redshift),1.0));
    printf("  Origin (x, y, z) = %lf, %lf, %lf\n", Origin_x, Origin_y, Origin_z);
    printf("  Nrep_min (x, y, z) = %d (%lf Mpc/h), %d (%lf Mpc/h), %d (%lf Mpc/h)\n", Nrep_neg_x, -Nrep_neg_x*Box-Origin_x, Nrep_neg_y, -Nrep_neg_y*Box-Origin_y, Nrep_neg_z, -Nrep_neg_z*Box-Origin_z);
    printf("  Nrep_max (x, y, z) = %d (%lf Mpc/h), %d (%lf Mpc/h), %d (%lf Mpc/h)\n", Nrep_pos_x, (Nrep_pos_x+1)*Box-Origin_x, Nrep_pos_y, (Nrep_pos_y+1)*Box-Origin_y, Nrep_pos_z, (Nrep_pos_z+1)*Box-Origin_z);
#endif
    printf("\nOutputs:\n");
    for (i=0; i<Noutputs; i++) printf("  Redshift = %lf, Nsteps = %d\n", OutputList[i].Redshift, OutputList[i].Nsteps);
    fflush(stdout);
  }   

  if (ThisTask == 0) {
    printf("\nInitialising Transfer Function/Power Spectrum\n");
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
  int j, k, m;
  int NoutputStart = 0;
  int timeStep, timeSteptot=0;
  unsigned int coord=0;
  double da=0;
  double A=1.0/(1.0+Init_Redshift);  // This is the scale factor which we'll be advancing below.
  double Di=growthD(A);              // initial growth factor
  double Di2=growthD2(A);            // initial 2nd order growth factor  
  double Dv=DprimeQ(A);              // T[D_{za}]=dD_{za}/dy
  double Dv2=growthD2v(A);           // T[D_{2lpt}]=dD_{2lpt}/dy

  // A is the scale factor of the particle positions.
  // AI is the scale factor of the particle velocities.
  // AF is the scale factor to which we should kick the particle velocities.
  // AFF is the scale factor to which we should drift the particle positions.
  double AI=A,AF=A,AFF=A;  

  displacement_fields();
    
  P = (struct part_data *) malloc((int)(ceil(NumPart*Buffer))*sizeof(struct part_data));

  // Generate the initial particle positions and velocities
  // If subtractLPT = 0 (non-COLA), then velocity is ds/dy, which is simply the 2LPT IC.
  // Else set vel = 0 if we subtract LPT. This is the same as the action of the operator L_- from TZE, as initial velocities are in 2LPT.
  for(i=0; i<Local_np; i++) {
    for (j=0; j<Nsample; j++) {
      for (k=0; k<Nsample; k++) {
        coord = (i * Nsample + j) * Nsample + k;
           
        P[coord].ID = ((unsigned long long)((i + Local_p_start) * Nsample + j)) * (unsigned long long)Nsample + (unsigned long long)k;
        for (m=0; m<3; m++) {
          P[coord].Dz[m] = ZA[m][coord];
          P[coord].D2[m] = LPT[m][coord];
          if (subtractLPT == 0) {
            P[coord].Vel[m] = P[coord].Dz[m]*Dv+P[coord].D2[m]*Dv2;
          } else {
            P[coord].Vel[m] = 0.0;
          }
        }

        P[coord].Pos[0] = periodic_wrap((i+Local_p_start)*(Box/(double)Nsample)+P[coord].Dz[0]*Di+P[coord].D2[0]*Di2);
        P[coord].Pos[1] = periodic_wrap(j*(Box/(double)Nsample)+P[coord].Dz[1]*Di+P[coord].D2[1]*Di2);
        P[coord].Pos[2] = periodic_wrap(k*(Box/(double)Nsample)+P[coord].Dz[2]*Di+P[coord].D2[2]*Di2);
      }
    }
  }

  for (i=0; i<3; i++) {
    free(ZA[i]);
    free(LPT[i]);
  }

  // If we want to output or start the lightcone at the initial redshift this is where we do it (it is tricky to compare
  // floating point numbers due to rounding errors so instead we see whether they are close)
  // ===================================================================================================================
  if (((Init_Redshift-OutputList[0].Redshift)/Init_Redshift <= 1.0E-6) || (Init_Redshift <= 1.0e-6)) {

#ifndef LIGHTCONE

    // Output particles.
    if (ThisTask == 0) {
      printf("Outputting Initial Conditions\n");
      printf("=============================\n\n");
    }

    sumx=0;
    sumy=0;
    sumz=0;

    Output(A,Dv,Dv2);

    // If this is the only output timestep then simply skip to the end
    if(Noutputs == 1) goto finalize;

#endif

    NoutputStart++;
  }

  // Now, we get to the N-Body part where we evolve with time via the Kick-Drift-Kick Method
  // =======================================================================================

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

  // Loop over all the timesteps in the timestep list
  // ================================================
#ifdef LIGHTCONE
  for (i=NoutputStart;i<Noutputs;i++) {
#else
  for (i=NoutputStart;i<=Noutputs;i++) {
#endif

    int nsteps=0;
    double ao=0;
    if (i == Noutputs) {
      nsteps = 1;
    } else {
      nsteps = OutputList[i].Nsteps;
      ao  = 1.0/(1.0+OutputList[i].Redshift);
      if (stepDistr == 0) da=(ao-A)/((double)nsteps);
      if (stepDistr == 1) da=(log(ao)-log(A))/((double)nsteps);
      if (stepDistr == 2) da=(CosmoTime(ao)-CosmoTime(A))/((double)nsteps);
    }

    // Perform the required number of timesteps between outputs
    // ========================================================
    for (timeStep=0;timeStep<nsteps;timeStep++) {

      timeSteptot++;
      if (ThisTask == 0) {
        printf("Iteration = %d\n------------------\n",timeSteptot);
        printf("a = %lf\n",A);
        printf("z = %lf\n",1.0/A-1.0);
        fflush(stdout);
      }

      // Calculate the particle accelerations for this timestep
      // ======================================================
      GetDisplacements();

      // Kick the particle velocities
      // ============================
      if (ThisTask == 0) {
        printf("Kicking the particles...\n");
        fflush(stdout);
      }

      /**********************************************************************************************/
      // If we wanted to interpolate the lightcone velocities we could put section currently in the // 
      // Drift subroutine here. This would allow us the have the velocities at AF and AFF which we  //
      // can use to invert the previous drift step and get the particle positions at AF and AFF     //                                                                                             //
      /**********************************************************************************************/

      // Half timestep for kick at output redshift to bring the velocities and positions to the same time
      if ((timeStep == 0) && (i != NoutputStart)) {
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

      Kick(AI,AF,A,Di);

#ifndef LIGHTCONE

      // If we are at an output timestep we modify the velocities and output the
      // particles. Then, if we are not yet at the end of the simulation, we update  
      // the velocity again up to the middle of the next timestep as per the usual KDK method.
      // =====================================================================================
      if ((timeStep == 0) && (i != NoutputStart)) {

        if (ThisTask == 0) {
          printf("Outputting the particles...\n");
          fflush(stdout);
        }

        // At the output timestep, add back LPT velocities if we had subtracted them. 
        // This corresponds to L_+ operator in TZE.
        Dv  = DprimeQ(A);    // dD_{za}/dy
        Dv2 = growthD2v(A);  // dD_{2lpt}/dy

        Output(A,Dv,Dv2);

        // If we have reached the last output timestep we skip to the end
        if(i == Noutputs) {
          if (ThisTask == 0) {
            printf("Iteration %d finished\n------------------\n\n", timeSteptot);
            fflush(stdout);
          }
          goto finalize;
        }

        // Otherwise we simply update the velocity to the middle of the timestep, where it
        // would have been if we weren't outputting. This involves only recalculating and applying `dda'
        // as the acceleration remains the same as calculated earlier
        AI = A;
        if (stepDistr == 0) AF=A+da*0.5;
        if (stepDistr == 1) AF=A*exp(da*0.5);
        if (stepDistr == 2) AF=AofTime((CosmoTime(AFF)+CosmoTime(A))*0.5); 

        sumDx=0;
        sumDy=0;
        sumDz=0;

        Kick(AI,AF,A,Di);      
      }

#endif

      for (j=0; j<3; j++) free(Disp[j]);

      // Drift the particle positions
      // ============================
      if (ThisTask == 0) {
        printf("Drifting the particles...\n");
        fflush(stdout);
      }

      if (stepDistr == 0) AFF=A+da;
      if (stepDistr == 1) AFF=A*exp(da);
      if (stepDistr == 2) AFF=AofTime(CosmoTime(A)+da);

#ifdef LIGHTCONE
      if (i == Noutputs-1) {
        Drift_Lightcone(A,AFF,AF,Di,timeStep);
      } else {
        Drift(A,AFF,AF,Di);
      }
#else
      Drift(A,AFF,AF,Di);
#endif


      // Step in time
      // ================
      A  = AFF;
      AI = AF; 

      Di = growthD(A);

      if (ThisTask == 0) {
        printf("Iteration %d finished\n------------------\n\n", timeSteptot);
        fflush(stdout);
      }
     
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Here is the last little bit
  // ===========================
#ifndef LIGHTCONE
  finalize:
#endif

  if (ThisTask == 0) {
    printf("Finishing up\n");
    printf("============\n");
    fflush(stdout);
  }

  free_powertable();
  free_transfertable();

  free(P);
  free(OutputList);
  free(Slab_to_task);
  free(Part_to_task);
  free(Local_nx_table);
  free(Local_np_table);
#ifdef GENERIC_FNL
  free(KernelTable);
#endif
#ifdef LIGHTCONE
  free(repflag);
#endif
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

// Kicking the particle velocities
// ===============================
void Kick(double AI, double AF, double A, double Di) {
  unsigned int n;
  double dda;
  double q1,q2;
  double ax,ay,az;
  double Om143=pow(Omega/(Omega+(1-Omega)*A*A*A),1./143.);
    
  sumx=0;
  sumy=0;
  sumz=0; 

  if (StdDA == 0) {
    dda=Sphi(AI,AF,A);
  } else if (StdDA == 1) {
    dda=(AF-AI)*A/Qfactor(A);
  } else {
    dda=SphiStd(AI,AF);
  }  
    
  q2=1.5*Omega*Di*Di*(1.0+7./3.*Om143)*A;      // T^2[D_{2lpt}]=d^2 D_{2lpt}/dy^2
  q1=1.5*Omega*Di*A;                           // T^2[D_{ZA}]=d^2 D_{ZA}/dy^2
    
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

  // Make sumx, sumy and sumz global averages
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  
    
  sumx /= (double)TotNumPart;  // We will subtract these to conserve momentum. 
  sumy /= (double)TotNumPart;  // Should be conserved, but just in case 3-linear interpolation makes a problem.
  sumz /= (double)TotNumPart;  // Never checked whether this makes a difference.
}

// Drifting the particle positions
// ===============================
void Drift(double A, double AFF, double AF, double Di) {

  unsigned int n;
  double dyyy;
  double da1,da2;
    
  if (StdDA == 0) {
    dyyy=Sq(A,AFF,AF);
  } else if (StdDA == 1) {
    dyyy=(AFF-A)/Qfactor(AF);
  } else {
    dyyy=SqStd(A,AFF);
  }

  da1=growthD(AFF)-Di;    // change in D
  da2=growthD2(AFF)-growthD2(A); // change in D_{2lpt}

  for(n=0; n<NumPart; n++) {
    P[n].Pos[0] += (P[n].Vel[0]-sumx)*dyyy;
    P[n].Pos[1] += (P[n].Vel[1]-sumy)*dyyy;
    P[n].Pos[2] += (P[n].Vel[2]-sumz)*dyyy;

    P[n].Pos[0] = periodic_wrap(P[n].Pos[0]+subtractLPT*(P[n].Dz[0]*da1+P[n].D2[0]*da2));
    P[n].Pos[1] = periodic_wrap(P[n].Pos[1]+subtractLPT*(P[n].Dz[1]*da1+P[n].D2[1]*da2));
    P[n].Pos[2] = periodic_wrap(P[n].Pos[2]+subtractLPT*(P[n].Dz[2]*da1+P[n].D2[2]*da2));
  }
}

// Output the data
// ===============
void Output(double A, double Dv, double Dv2) {

  FILE * fp; 
  char buf[300];
  int nprocgroup, groupTask, masterTask;
  unsigned int n;
  double Z = (1.0/A)-1.0;
  double fac = Hubble/pow(A,1.5);
  double lengthfac = UnitLength_in_cm/3.085678e24;     // Convert positions to Mpc/h
  double velfac    = UnitVelocity_in_cm_per_s/1.0e5;   // Convert velocities to km/s

  // Remember to add the ZA and 2LPT velocities back on and convert to PTHalos velocity units
#ifdef GADGET_STYLE
  size_t bytes;
  int k, pc, dummy, blockmaxlen;
  float * block;

  for (n=0; n<NumPart; n++) {
    P[n].Pos[0] *= lengthfac;
    P[n].Pos[1] *= lengthfac;
    P[n].Pos[2] *= lengthfac;

    P[n].Vel[0] = velfac*fac*(P[n].Vel[0]-sumx+(P[n].Dz[0]*Dv+P[n].D2[0]*Dv2)*subtractLPT);
    P[n].Vel[1] = velfac*fac*(P[n].Vel[1]-sumy+(P[n].Dz[1]*Dv+P[n].D2[1]*Dv2)*subtractLPT);
    P[n].Vel[2] = velfac*fac*(P[n].Vel[2]-sumz+(P[n].Dz[2]*Dv+P[n].D2[2]*Dv2)*subtractLPT);
  }
#endif

  nprocgroup = NTask / NumFilesWrittenInParallel;
  if (NTask % NumFilesWrittenInParallel) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if(NumPart > 0) {
        sprintf(buf, "%s/%s_z%dp%03d.%d", OutputDir, FileBase, (int)Z, (int)rint((Z-(int)Z)*1000), ThisTask);
        if(!(fp = fopen(buf, "w"))) {
          printf("\nERROR: Can't write in file '%s'.\n\n", buf);
          FatalError("main.c", 597);
        }
#ifdef GADGET_STYLE
        // Gadget header stuff
        for(k = 0; k < 6; k++) {
          header.npart[k] = 0;
          header.npartTotal[k] = 0;
          header.mass[k] = 0;
        }
        header.npart[1] = NumPart;
        header.npartTotal[1] = TotNumPart;
        header.npartTotal[2] = (TotNumPart >> 32);
        header.mass[1] = (3.0*Omega*Hubble*Hubble*Box*Box*Box) / (8.0*PI*G*TotNumPart);
        header.time = A;
        header.redshift = Z;

        header.flag_sfr = 0;
        header.flag_feedback = 0;
        header.flag_cooling = 0;
        header.flag_stellarage = 0;
        header.flag_metals = 0;
        header.flag_stellarage = 0;
        header.flag_metals = 0;
        header.hashtabsize = 0;

        header.num_files = NTaskWithN;

        header.BoxSize = Box;
        header.Omega0 = Omega;
        header.OmegaLambda = 1.0-Omega;
        header.HubbleParam = HubbleParam;

        dummy = sizeof(header);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        my_fwrite(&header, sizeof(header), 1, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        block = (float *)malloc(bytes = 10 * 1024 * 1024);
        blockmaxlen = bytes / (3 * sizeof(float));

        // write coordinates
        dummy = sizeof(float) * 3 * NumPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          for(k = 0; k < 3; k++) block[3 * pc + k] = (float)P[n].Pos[k];
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3 * pc, fp);
	    pc = 0;
	  }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3 * pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        // write velocities
        dummy = sizeof(float) * 3 * NumPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          for(k = 0; k < 3; k++) block[3 * pc + k] = (float)P[n].Vel[k];
          pc++;
          if(pc == blockmaxlen) {
	    my_fwrite(block, sizeof(float), 3 * pc, fp);
	    pc = 0;
	  }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3 * pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        free(block);   
#else
        for(n=0; n<NumPart; n++){
          double P_Vel[3];
          P_Vel[0] = fac*(P[n].Vel[0]-sumx+(P[n].Dz[0]*Dv+P[n].D2[0]*Dv2)*subtractLPT);
          P_Vel[1] = fac*(P[n].Vel[1]-sumy+(P[n].Dz[1]*Dv+P[n].D2[1]*Dv2)*subtractLPT);
          P_Vel[2] = fac*(P[n].Vel[2]-sumz+(P[n].Dz[2]*Dv+P[n].D2[2]*Dv2)*subtractLPT);
          fprintf(fp,"%12llu %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                      P[n].ID, (float)(lengthfac*P[n].Pos[0]),(float)(lengthfac*P[n].Pos[1]),(float)(lengthfac*P[n].Pos[2]),(float)(velfac*P_Vel[0]),(float)(velfac*P_Vel[1]),(float)(velfac*P_Vel[2]));
        }
#endif
        fclose(fp);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

#ifdef GADGET_STYLE
  for (n=0; n<NumPart; n++) {
    P[n].Pos[0] /= lengthfac;
    P[n].Pos[1] /= lengthfac;
    P[n].Pos[2] /= lengthfac;

    P[n].Vel[0] = P[n].Vel[0]/(velfac*fac)+sumx-(P[n].Dz[0]*Dv+P[n].D2[0]*Dv2)*subtractLPT;
    P[n].Vel[1] = P[n].Vel[1]/(velfac*fac)+sumy-(P[n].Dz[1]*Dv+P[n].D2[1]*Dv2)*subtractLPT;
    P[n].Vel[2] = P[n].Vel[2]/(velfac*fac)+sumz-(P[n].Dz[2]*Dv+P[n].D2[2]*Dv2)*subtractLPT;
  }
#endif  
 
  return;
}
