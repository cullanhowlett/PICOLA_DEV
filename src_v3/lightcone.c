
#include "vars.h"
#include "proto.h"

// Set up the lightcone parameters
// ===============================
void set_lightcone(void) {

  Light = LIGHT / UnitLength_in_cm * UnitTime_in_s;

  // Calculate the actual maximum number of replicate boxes we need in all six directions. This allows
  // us to remove replicates that are unnecessary
  double Rcomov_max = Light/Hubble*SphiStd(1.0/(1.0+OutputList[0].Redshift),1.0);
  if ((int)(ceil(fabs((Origin_x - Rcomov_max)/Box))) < Nrep_neg_x) Nrep_neg_x = (int)(ceil(fabs((Origin_x - Rcomov_max)/Box)));
  if ((int)(ceil(fabs((Origin_y - Rcomov_max)/Box))) < Nrep_neg_y) Nrep_neg_y = (int)(ceil(fabs((Origin_y - Rcomov_max)/Box)));
  if ((int)(ceil(fabs((Origin_z - Rcomov_max)/Box))) < Nrep_neg_z) Nrep_neg_z = (int)(ceil(fabs((Origin_z - Rcomov_max)/Box)));
  if ((int)(ceil((Rcomov_max + Origin_x)/Box))-1 < Nrep_pos_x) Nrep_pos_x = (int)(ceil((Rcomov_max + Origin_x)/Box))-1;
  if ((int)(ceil((Rcomov_max + Origin_y)/Box))-1 < Nrep_pos_y) Nrep_pos_y = (int)(ceil((Rcomov_max + Origin_y)/Box))-1;
  if ((int)(ceil((Rcomov_max + Origin_z)/Box))-1 < Nrep_pos_z) Nrep_pos_z = (int)(ceil((Rcomov_max + Origin_z)/Box))-1;

  // Create an array of flags that for each replicate tell us whether or not we need to check the particles in it.
  // A flag of 0 means we check it, 1 means the box is completely inside the lightcone so no need to check it this step 
  // 2 means the box has completely left the lightcone so no need to ever check it again.
  repflag = (int *)calloc((Nrep_neg_x+Nrep_pos_x+1)*(Nrep_neg_y+Nrep_pos_y+1)*(Nrep_neg_z+Nrep_pos_z+1),sizeof(int));
  
  Nrep_neg_max[0] = Nrep_neg_x; Nrep_pos_max[0] = Nrep_pos_x;
  Nrep_neg_max[1] = Nrep_neg_y; Nrep_pos_max[1] = Nrep_pos_y;
  Nrep_neg_max[2] = Nrep_neg_z; Nrep_pos_max[2] = Nrep_pos_z;
  

  return;
}

// Drift and output the particles for lightcone simulations
// ========================================================
void Drift_Lightcone(double A, double AFF, double AF, double Di, int timeStep) {

  // We'll flag to see if the particle has left the lightcone here and output if it has
  // We don't bother interpolating the velocity, only the position, as the implicit assumption
  // with KDK anyway is that the velocity at the halfway point in the timestep is constant
  // between the initial and final particle positions for that timestep. If we did want to interpolate
  // the velocity, we could move this whole section to the location marked above, then reverse the 
  // updating and periodic wrapping of the particle positions to allow interpolation.
  // To avoid having to store all the replicate of the particles, we also output the particles as soon 
  // as they are flagged, allowing us to loop over all the replicated boxes instead. Finally, because all
  // the lightcone stuff is done here, we don't need memory to store the particle flags.

  FILE * fp; 
  char buf[300];
  int i, j, k, ii, jj, kk;
  int flag, nprocgroup, groupTask, masterTask;
  int repcount_low, repcount_high, coord;
  unsigned int n;
  double dyyy, da1, da2, dv1, dv2;
  double dyyy_tmp, da1_tmp, da2_tmp;
  double Delta_Pos[3], P_Pos[3], P_Vel[3];
  double boundary = 10.0;
  double fac = Hubble/pow(AF,1.5);
  double lengthfac = UnitLength_in_cm/3.085678e24;     // Convert positions to Mpc/h
  double velfac    = UnitVelocity_in_cm_per_s/1.0e5;   // Convert velocities to km/s
  double Rcomov_old  = Light/Hubble*SphiStd(A,1.0);
  double Rcomov_new  = Light/Hubble*SphiStd(AFF,1.0);
  double Rcomov_old2 = Rcomov_old*Rcomov_old; 
  double Rcomov_new2 = Rcomov_new*Rcomov_new;
  double Xpart, Ypart, Zpart, Rpart_old, Rpart_new, AL;   

  if (StdDA == 0) {
    dyyy=Sq(A,AFF,AF);
  } else if (StdDA == 1) {
    dyyy=(AFF-A)/Qfactor(AF);
  } else {
    dyyy=SqStd(A,AFF);
  }

  da1=growthD(AFF)-Di;           // change in D
  da2=growthD2(AFF)-growthD2(A); // change in D_{2lpt}

  // Add back LPT velocities if we had subtracted them. 
  // This corresponds to L_+ operator in TZE.
  dv1 = DprimeQ(AF);    // dD_{za}/dy
  dv2 = growthD2v(AF);  // dD_{2lpt}/dy

  // Check the maximum number of replicates in each direction and update if necessary
  if ((int)(ceil(fabs((Origin_x - Rcomov_old)/Box))) < Nrep_neg_x) Nrep_neg_x = (int)(ceil(fabs((Origin_x - Rcomov_old)/Box)));
  if ((int)(ceil(fabs((Origin_y - Rcomov_old)/Box))) < Nrep_neg_y) Nrep_neg_y = (int)(ceil(fabs((Origin_y - Rcomov_old)/Box)));
  if ((int)(ceil(fabs((Origin_z - Rcomov_old)/Box))) < Nrep_neg_z) Nrep_neg_z = (int)(ceil(fabs((Origin_z - Rcomov_old)/Box)));
  if ((int)(ceil((Rcomov_old + Origin_x)/Box))-1 < Nrep_pos_x) Nrep_pos_x = (int)(ceil((Rcomov_old + Origin_x)/Box))-1;
  if ((int)(ceil((Rcomov_old + Origin_y)/Box))-1 < Nrep_pos_y) Nrep_pos_y = (int)(ceil((Rcomov_old + Origin_y)/Box))-1;
  if ((int)(ceil((Rcomov_old + Origin_z)/Box))-1 < Nrep_pos_z) Nrep_pos_z = (int)(ceil((Rcomov_old + Origin_z)/Box))-1;

  if (ThisTask == 0) printf("%lf, %lf, %d, %d, %d, %d, %d, %d\n", Rcomov_old, Rcomov_new, Nrep_neg_x, Nrep_pos_x, Nrep_neg_y, Nrep_pos_y, Nrep_neg_z, Nrep_pos_z);

  // Loop over all replicates and if any replicate is completely within/outside the lightcone (i.e. all eight vertices are within Rcomov_new/outside Rcomov_old)
  // then flag it so that we don't have to loop over it. NOTE: this method will not work to see if replicates are completely outside the
  // lightcone if we also introduce an angular restriction on the lightcone, as I then can easily think of a way to place a lightcone so that 
  // all eight vertices of a replicated box are outside the lightcone, yet the lightcone still passes through it. In this case the only way 
  // I can think of to check if a box is completely outside the lightcone is to loop over all the particles once and flag it if no particles 
  // are within Rcomov_old.
  for (i = -Nrep_neg_x; i<=Nrep_pos_x; i++) {
    for (j = -Nrep_neg_y; j<=Nrep_pos_y; j++) {
      for (k = -Nrep_neg_z; k<=Nrep_pos_z; k++) {

        coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);

        // Skip this replicate if we already know it is completely outside the lightcone
        // the particle positions first
        if (repflag[coord] == 2) continue;
 
        // Loop over all the vertices
        repflag[coord] = 0;
        repcount_low = repcount_high = 0;
        for (ii = 0; ii < 2; ii++) {
          for (jj = 0; jj < 2; jj++) {
            for (kk = 0; kk < 2; kk++) {

              Xpart = i*Box - Origin_x + (ThisTask+ii)*(Box/(double)NTask);
              Ypart = (j+jj)*Box - Origin_y;
              Zpart = (k+kk)*Box - Origin_z;
              Rpart_old = Xpart*Xpart+Ypart*Ypart+Zpart*Zpart;
                  
              // Include a buffer region to account for the fact that the particle might move beyond the box boundaries, 500kpc should be enough.
              if (Rpart_old <= Rcomov_new2-boundary) {
                repcount_low++;
              } else if (Rpart_old >= Rcomov_old2+boundary) {
                repcount_high++;
              }
            }
          }
        }
            
        // If box is completely inside or outside the lightcone, we skip it UNLESS it is the last replicate, in which case update
        // the particle positions first
        if (repcount_low == 8) {
          repflag[coord] = 1;
        } else if (repcount_high == 8) {
          repflag[coord] = 2;
        }
      }
    }
  }

  // Loop over all the processors that are allowed to output at once
  nprocgroup = NTask / NumFilesWrittenInParallel;
  if (NTask % NumFilesWrittenInParallel) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if(NumPart > 0) {
        sprintf(buf, "%s/%s_lightcone.%d", OutputDir, FileBase, ThisTask);
        if (timeStep == 0) {
          // Overwrite any pre-existing output files otherwise we'll append onto the end of them.
          if(!(fp = fopen(buf, "w"))) {
            printf("\nERROR: Can't write in file '%s'.\n\n", buf);
            FatalError("lightcone.c", 93);
          }
        } else {
          if(!(fp = fopen(buf, "a"))) {
            printf("\nERROR: Can't write in file '%s'.\n\n", buf);
            FatalError("lightcone.c", 98);
          }
        }

        // Loop over all particles, modifying the position based on the current replicate
        for(n=0; n<NumPart; n++) {

          if (n % 100000 == 0) printf("%d, %u\n", ThisTask, n);

          Delta_Pos[0] = (P[n].Vel[0]-sumx)*dyyy+subtractLPT*(P[n].Dz[0]*da1+P[n].D2[0]*da2);
          Delta_Pos[1] = (P[n].Vel[1]-sumy)*dyyy+subtractLPT*(P[n].Dz[1]*da1+P[n].D2[1]*da2);   
          Delta_Pos[2] = (P[n].Vel[2]-sumz)*dyyy+subtractLPT*(P[n].Dz[2]*da1+P[n].D2[2]*da2);     

          // Check that 500kpc boundaries is enough
          if (Delta_Pos[0]*Delta_Pos[0]+Delta_Pos[1]*Delta_Pos[1]+Delta_Pos[2]*Delta_Pos[2] > boundary) {
            printf("\nERROR: Particle displacement greater than boundary for lightcone replicate estimate.\n");
            printf("       increase boundary condition in lightcone.c (line 56)\n\n");
            FatalError("lightcone.c", 212);
          }

          // Loop over all replicates
          for (i = -Nrep_neg_x; i<=Nrep_pos_x; i++) {
            for (j = -Nrep_neg_y; j<=Nrep_pos_y; j++) {
              for (k = -Nrep_neg_z; k<=Nrep_pos_z; k++) {

                coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);
                if (repflag[coord] != 0) continue;

                // Did the particle start the timestep inside the lightcone?
                flag = 0;
                Xpart = P[n].Pos[0] - Origin_x + (i*Box);
                Ypart = P[n].Pos[1] - Origin_y + (j*Box);
                Zpart = P[n].Pos[2] - Origin_z + (k*Box);
                Rpart_old = Xpart*Xpart+Ypart*Ypart+Zpart*Zpart;
 
                if (Rpart_old <= Rcomov_old2) flag = 1;

                // Have any particles that started inside the lightcone now exited?
                if (flag) {
                  Xpart = P[n].Pos[0] + Delta_Pos[0] - Origin_x + (i*Box);
                  Ypart = P[n].Pos[1] + Delta_Pos[1] - Origin_y + (j*Box);
                  Zpart = P[n].Pos[2] + Delta_Pos[2] - Origin_z + (k*Box);
                  Rpart_new = Xpart*Xpart+Ypart*Ypart+Zpart*Zpart;

                  if (Rpart_new > Rcomov_new2) {
  
                    // Interpolate the particle position. We do this by first calculating the exact time at which
                    // the particle exited the lightcone, then updating the position to there.
                    AL = A + (AFF-A)*((Rcomov_old-sqrt(Rpart_old))/((sqrt(Rpart_new)-sqrt(Rpart_old))-(Rcomov_new-Rcomov_old)));

                    if (StdDA == 0) {
                      dyyy_tmp=Sq(A,AL,AF);
                    } else if (StdDA == 1) {
                      dyyy_tmp=(AL-A)/Qfactor(AF);
                    } else {
                      dyyy_tmp=SqStd(A,AL);
                    }

                    da1_tmp=growthD(AL)-Di;           // change in D
                    da2_tmp=growthD2(AL)-growthD2(A); // change in D_{2lpt}
              
                    P_Pos[0] = P[n].Pos[0] + (P[n].Vel[0]-sumx)*dyyy_tmp+subtractLPT*(P[n].Dz[0]*da1_tmp+P[n].D2[0]*da2_tmp) + (i*Box);
                    P_Pos[1] = P[n].Pos[1] + (P[n].Vel[1]-sumy)*dyyy_tmp+subtractLPT*(P[n].Dz[1]*da1_tmp+P[n].D2[1]*da2_tmp) + (j*Box);
                    P_Pos[2] = P[n].Pos[2] + (P[n].Vel[2]-sumz)*dyyy_tmp+subtractLPT*(P[n].Dz[2]*da1_tmp+P[n].D2[2]*da2_tmp) + (k*Box);

                    P_Vel[0] = fac*(P[n].Vel[0]-sumx+(P[n].Dz[0]*dv1+P[n].D2[0]*dv2)*subtractLPT);
                    P_Vel[1] = fac*(P[n].Vel[1]-sumy+(P[n].Dz[1]*dv1+P[n].D2[1]*dv2)*subtractLPT);
                    P_Vel[2] = fac*(P[n].Vel[2]-sumz+(P[n].Dz[2]*dv1+P[n].D2[2]*dv2)*subtractLPT);

                    fprintf(fp,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                                (float)(lengthfac*P_Pos[0]),(float)(lengthfac*P_Pos[1]),(float)(lengthfac*P_Pos[2]),(float)(velfac*P_Vel[0]),(float)(velfac*P_Vel[1]),(float)(velfac*P_Vel[2]));
                  }
                }
              }
            }
          }
 
          // Update the particle's position
          P[n].Pos[0] = periodic_wrap(P[n].Pos[0]+Delta_Pos[0]);
          P[n].Pos[1] = periodic_wrap(P[n].Pos[1]+Delta_Pos[1]);
          P[n].Pos[2] = periodic_wrap(P[n].Pos[2]+Delta_Pos[2]); 
        }
        fclose(fp);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
 
  return;
}
