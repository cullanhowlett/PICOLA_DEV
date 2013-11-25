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

/* ==========================================================================*/
/* v1.0: This file contains all the prototypes for function used in the code.*/
/* ==========================================================================*/

// main.c
void Output(double A, double Dv, double Dv2);
void Kick(double AI, double AF, double A, double Di);
void Drift(double A, double AFF, double AF, double Di);

// cosmo.c
double gpQ(double a);
double DERgpQ(double a);
double decayD(double a);
double DprimeQ(double a);
double Qfactor(double a);
double AofTime(double y);
double growthD(double a);
double growthD2(double a);
double growthD2v(double a);
double CosmoTime(double af);
double growthDtemp(double a);
double growthD2temp(double a);
double SqStd(double ai,double af);
double SphiStd(double ai,double af);
double fun(double x, void * params);
double funSqStd(double a, void * params);
double AofTimeFun(double a,void * params);
double Sq(double ai,double af,double aRef);
double funSphiStd(double a, void * params);
double Sphi(double ai,double af,double aRef);
double CosmoTimeFun (double a, void * params);
double GrowthFactor(double astart, double aend);

// auxPM.c
void Forces(void);
void PtoMesh(void);
void MtoParticles(void);
void MoveParticles(void);
void GetDisplacements(void);
void FatalError(char * filename, int linenum);
#if (MEMORY_MODE || SINGLE_PRECISION)
float periodic_wrap(float x);
#else
double periodic_wrap(double x);
#endif
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);

// read_param.c
void read_outputs(void);
void read_parameterfile(char * fname);
int sort_redshift(const void * Item1, const void * Item2);

// kernel.c
#ifdef GENERIC_FNL
void read_kernel_table(void);
#endif

// lightcone.c
#ifdef LIGHTCONE
void set_lightcone(void);
void Drift_Lightcone(double A, double AFF, double AF, double Di, int timeStep);
#endif

// 2LPT.c
void set_units(void);
void initialize_ffts(void);
void initialize_parts(void);
void displacement_fields(void);

// power.c
void print_spec(void);
void free_powertable(void);
void read_power_table(void);
void free_transfertable(void);
void read_transfer_table(void);
void initialize_powerspectrum(void);
void initialize_transferfunction(void);
int compare_logk(const void *a, const void *b);
int compare_transfer_logk(const void *a, const void *b);
double fnl(double x);
double tk_eh(double k);
double TransferFunc(double k);
double PowerSpec(double kmag);
double TopHatSigma2(double R);
double PowerSpec_EH(double k);
double tk_Efstathiou(double k);
double TransferFunc_EH(double k);
double PowerSpec_Tabulated(double k);
double PowerSpec_Efstathiou(double k);
double TransferFunc_Tabulated(double k);
double TransferFunc_Efstathiou(double k);
double sigma2_int(double k, void * params);
