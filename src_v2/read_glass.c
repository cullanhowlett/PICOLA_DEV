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

/* ======================================================================================================*/
/* v1.0: This file contains the routines required to create the initial glass-like particle distribution.*/
/* v2.0: Unchanged.                                                                                      */
/* ======================================================================================================*/

#include "vars.h"
#include "proto.h"

void read_glass(char * fname) {

  FILE *fd = 0;
  char buf[500];
  int num,numfiles,skip;
  int i,j,k,n,slab,type;
  unsigned int *npart_Task;
  unsigned int m, count, nlocal, dummy, dummy2;
  unsigned long long IDStart;
  float *pos = 0;
  float x, y, z;
  double buffer = 1.1;          // This is the amount of extra space to include for particles that might move between tasks during timesteps
  size_t bytes;

#define SKIP {my_fread(&dummy, sizeof(int), 1, fd);}
#define SKIP2 {my_fread(&dummy2, sizeof(int), 1, fd);}

  if(ThisTask == 0) {
    printf("Reading Lagrangian glass file...\n");
    fflush(stdout);

    numfiles = find_files(fname);

    for(num = 0, skip = 0; num < numfiles; num++) {

      if(numfiles > 1) {
        sprintf(buf, "%s.%d", fname, num);
      } else {
	sprintf(buf, "%s", fname);
      }

      if(!(fd = fopen(buf, "r"))) {
        printf("\nERROR: Can't open glass file '%s'.\n", buf);
	FatalError(1);
      }

      SKIP;
      my_fread(&header1, sizeof(header1), 1, fd);
      SKIP2;

      if(dummy != sizeof(header1) || dummy2 != sizeof(header1)) {
        printf("\nERROR: Incorrect header size.\n");
	FatalError(2);
      }

      nlocal = 0;
      for(k = 0; k < 6; k++) nlocal += header1.npart[k];

      printf("Reading '%s' with %d particles...\n", fname, nlocal);

      if(num == 0) {
        Nglass = 0;
	for(k = 0; k < 6; k++) Nglass += header1.npartTotal[k];

	pos = (float *) malloc(sizeof(float) * Nglass * 3);

	if(!(pos)) {
          printf("\nERROR: Failed to allocate %f MB on Task %d for glass file.\n",sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0),ThisTask);
          FatalError(112);
	}
      }

      SKIP;
      my_fread(&pos[3 * skip], sizeof(float), 3 * nlocal, fd);
      SKIP2;
 
      if(dummy != sizeof(float) * 3 * nlocal || dummy2 != sizeof(float) * 3 * nlocal) {
        printf("\nERROR: Incorrect block structure in positions block.\n");
	FatalError(3);
      }
      skip += nlocal;

      fclose(fd);
    }
  }

  MPI_Bcast(&Nglass, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&header1, sizeof(header1), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask != 0) {
    pos = (float *) malloc(sizeof(float) * Nglass * 3);

    if(!(pos)) {
      printf("\nERROR: Failed to allocate %f MB on Task %d for glass file\n",sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0),ThisTask);
      FatalError(112);
    }
  }

  MPI_Bcast(&pos[0], sizeof(float) * Nglass * 3, MPI_BYTE, 0, MPI_COMM_WORLD);

  npart_Task = (unsigned int *)malloc(sizeof(unsigned int) * NTask);
  for(i = 0; i < NTask; i++) npart_Task[i] = 0;

  for(i = 0; i < GlassTileFac; i++) {
    for(j = 0; j < GlassTileFac; j++) {
      for(k = 0; k < GlassTileFac; k++) {
        for(type = 0, n = 0; type < 6; type++) {
	  for(m = 0; m < header1.npartTotal[type]; m++, n++) {

            x = pos[3 * n] / header1.BoxSize * (Box / GlassTileFac) + i * (Box / GlassTileFac);

	    slab = (int)(x / Box * Nmesh);
	    if(slab >= Nmesh) slab = Nmesh - 1;

	    npart_Task[Slab_to_task[slab]] += 1;
	  }
	}
      }
    }
  }

  TotNumPart = 0;
  NTaskWithN = 0;		
  NumPart = npart_Task[ThisTask];
  for(i = 0; i < NTask; i++) {
    TotNumPart += (unsigned long long)npart_Task[i];
    if(npart_Task[i] > 0) NTaskWithN++;
  }

  if(ThisTask == 0) {
    printf("\nParticles\n---------------------\n");
    for(i = 0; i < NTask; i++) printf("Task = %d: Particles = %d\n", i, npart_Task[i]);
    printf("\n----------------------\n");
    printf("Total number of particles = %llu\n\n", TotNumPart);
    fflush(stdout);
  }
  free(npart_Task);

  NumPartWithBuf = (unsigned int)rint(NumPart * buffer);
  if(NumPart) {
    P = (struct part_data *) malloc(bytes = sizeof(struct part_data) * NumPartWithBuf);
    if(!(P)) {
      printf("\nERROR: Failed to allocate %g MB (%d particles) on Task %d\n",bytes / (1024.0 * 1024.0),NumPartWithBuf,ThisTask);
      FatalError(9891);
    }
  }

  count = 0;
  IDStart = 1;
  for(i = 0; i < GlassTileFac; i++) {
    for(j = 0; j < GlassTileFac; j++) {
      for(k = 0; k < GlassTileFac; k++) {
        for(type = 0, n = 0; type < 6; type++) {
	  for(m = 0; m < header1.npartTotal[type]; m++, n++) {
	
            x = pos[3 * n] / header1.BoxSize * (Box / GlassTileFac) + i * (Box / GlassTileFac);

	    slab = (int)(x / Box * Nmesh);
	    if(slab >= Nmesh) slab = Nmesh - 1;

            if(Slab_to_task[slab] == ThisTask) {
	      y = pos[3 * n + 1] / header1.BoxSize * (Box / GlassTileFac) + j * (Box / GlassTileFac);
	      z = pos[3 * n + 2] / header1.BoxSize * (Box / GlassTileFac) + k * (Box / GlassTileFac);

	      P[count].Pos[0] = (float_kind)x;
	      P[count].Pos[1] = (float_kind)y;
	      P[count].Pos[2] = (float_kind)z;
	      P[count].ID = IDStart;

	      count++;
	    }
            IDStart++;
	  }
	}
      }
    }
  }

  if(count != NumPart) {
    printf("\nERROR: Fatal mismatch (%u %u) on Task %d\n", count, NumPart, ThisTask);
    FatalError(1);
  }
  free(pos);

  return;
}

int find_files(char * fname) {
  
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd = fopen(buf, "r"))) {
    my_fread(&dummy, sizeof(dummy), 1, fd);
    my_fread(&header, sizeof(header), 1, fd);
    my_fread(&dummy, sizeof(dummy), 1, fd);
    fclose(fd);
    return header.num_files;
  }

  if((fd = fopen(buf1, "r"))) {
    my_fread(&dummy, sizeof(dummy), 1, fd);
    my_fread(&header, sizeof(header), 1, fd);
    my_fread(&dummy, sizeof(dummy), 1, fd);
    fclose(fd);
    header.num_files = 1;
    return header.num_files;
  }

  FatalError(121);
  return 0;
}
