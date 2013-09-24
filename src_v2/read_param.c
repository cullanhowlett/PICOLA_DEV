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

/* =====================================================================================*/
/* v1.0: This file contains routines to read in the run parameters file.                */
/* v2.0: Modified to read in extra parameters for non-gaussianity and transfer function.*/
/*       Also now runs some checks to ensure compatible run parameters.                 */
/* =====================================================================================*/

#include "vars.h"
#include "proto.h"

void read_parameterfile(char * fname) {

#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  char buf[200],buf1[200],buf2[200],buf3[200];
  int i,j,nt;
  int id[MAXTAGS];
  int errorFlag = 0;

  // read parameter file on all processes for simplicity

  nt = 0;

  strcpy(tag[nt], "UseCOLA");
  addr[nt] = &UseCOLA;
  id[nt++] = INT;

  strcpy(tag[nt], "nsteps");
  addr[nt] = &nsteps;
  id[nt++] = INT;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaBaryon");
  addr[nt] = &OmegaBaryon;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "HubbleParam");
  addr[nt] = &HubbleParam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ShapeGamma");
  addr[nt] = &ShapeGamma;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "PrimordialIndex");
  addr[nt] = &PrimordialIndex;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Box");
  addr[nt] = &Box;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Buffer");
  addr[nt] = &Buffer;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Init_Redshift");
  addr[nt] = &Init_Redshift;
  id[nt++] = FLOAT;

#ifndef GAUSSIAN
  strcpy(tag[nt], "Fnl_Redshift");
  addr[nt] = &Fnl_Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Fnl");
  addr[nt] = &Fnl;
  id[nt++] = FLOAT;
#endif

  strcpy(tag[nt], "Final_Redshift");
  addr[nt] = &Final_Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nmesh");
  addr[nt] = &Nmesh;
  id[nt++] = INT;

  strcpy(tag[nt], "Nsample");
  addr[nt] = &Nsample;
  id[nt++] = INT;

  strcpy(tag[nt], "FileWithInputSpectrum");
  addr[nt] = FileWithInputSpectrum;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputTransfer");
  addr[nt] = FileWithInputTransfer;
  id[nt++] = STRING;

#ifdef GENERIC_FNL
  strcpy(tag[nt], "FileWithInputKernel");
  addr[nt] = FileWithInputKernel;
  id[nt++] = STRING;
#endif

  strcpy(tag[nt], "Seed");
  addr[nt] = &Seed;
  id[nt++] = INT;

  strcpy(tag[nt], "SphereMode");
  addr[nt] = &SphereMode;
  id[nt++] = INT;

  strcpy(tag[nt], "NumFilesWrittenInParallel");
  addr[nt] = &NumFilesWrittenInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileBase");
  addr[nt] = FileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "WhichSpectrum");
  addr[nt] = &WhichSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "WhichTransfer");
  addr[nt] = &WhichTransfer;
  id[nt++] = INT;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
  addr[nt] = &InputSpectrum_UnitLength_in_cm;
  id[nt++] = FLOAT;

  if((fd = fopen(fname, "r"))) {
    while(!feof(fd)) {
      buf[0] = 0;
      fgets(buf, 200, fd);

      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
      if(buf1[0] == '%') continue;

      for(i = 0, j = -1; i < nt; i++) {
        if(strcmp(buf1, tag[i]) == 0)  {
          j = i;
	  tag[i][0] = 0;
	  break;
	}
      }
      
      if(j >= 0) {
	switch (id[j]) {
	  case FLOAT:
	    *((double *) addr[j]) = atof(buf2);
	    break;
	  case STRING:
	    strcpy((char *)addr[j], buf2);
	    break;
	  case INT:
	    *((int *) addr[j]) = atoi(buf2);
	    break;
	}
      } else {
        if(ThisTask == 0) fprintf(stdout,"\nERROR: In file %s:  Tag '%s' not allowed or multiple defined.\n",fname,buf1);
	errorFlag = 1;
      }
    }
    fclose(fd);
    for(i = 0; i < nt; i++) {
      if(*tag[i]) {
        if(ThisTask == 0) fprintf(stdout, "\nERROR: I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
        errorFlag = 1;
      }
    }
  } else {
    if(ThisTask == 0) fprintf(stdout,"\nERROR: Parameter file '%s' not found.\n",fname);
    errorFlag = 1;
  }

  // Check the run parameters to ensure compatible gaussian/non-gaussian options
  if((WhichSpectrum != 0) && (WhichTransfer !=0)) {
    if (ThisTask == 0) {
      printf("\nERROR: You are running with both power spectrum and tranfer function.\n");
      printf("         Please select the appropriate one.\n");
    }
    errorFlag = 1;
  }

#ifndef GAUSSIAN 
  if((WhichSpectrum != 0) || (WhichTransfer == 0)) { 
    if (ThisTask == 0) {
      printf("\nERROR: Non-Gaussian models require the transfer function as input.\n");
      printf("         Switch WhichSpectrum to zero in the input parameter file.\n");
    } 
    errorFlag = 1;
  }
#endif
#ifdef LOCAL_FNL 
   if(PrimordialIndex != 1.0) {
     if (ThisTask == 0) printf("\nERROR: Local non-gaussianity with tilted power spectrum requires the GENERIC_FNL option in the Makefile\n");
     errorFlag = 1;
   }
#endif 
#ifdef ORTHO_FNL 
   if(PrimordialIndex != 1.0) {
     if (ThisTask == 0) printf("\nERROR: Orthogonal non-gaussianity with tilted power spectrum requires the GENERIC_FNL option in the Makefile\n"); 
     errorFlag = 1;
   }
#endif 
#ifdef EQUIL_FNL 
   if(PrimordialIndex != 1.0) {
     if (ThisTask == 0) printf("\nERROR: Equilateral non-gaussianity with tilted power spectrum requires the GENERIC_FNL option in the Makefile\n");
     errorFlag = 1;
   }
#endif 

  if(errorFlag) {
    MPI_Finalize();
    exit(0);
  }

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS

  return;
}
