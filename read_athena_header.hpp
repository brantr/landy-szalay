#ifndef READ_ATHENA_HEADER_H
#define READ_ATHENA_HEADER_H
struct AthenaHeaderDef
{
#ifdef ATHENA4
  int CoordinateSystem;
#endif /* ATHENA4 */
  int nx;
  int ny;
  int nz;
  
  //; nvar=4 means isothermal hydro.  nvar=5 means adiabatic hydro
  //; nvar=7 means isothermal MHD.    nvar=8 means adiabatic mhd
  int nvar;
  int nscalars;
  int ngrav;

#ifdef ATHENA4
  int flag_tracers;
#endif /*ATHENA4*/
  /*int flag_gravity;
  int flag_scalars;
  int flag_mhd;
  int flag_isothermal;
  int flag_pressure;
	
  int ipressure;*/

  //(gamma-1) and isothermal sound speed

  float gamma_minus_1;	
  float c_s_iso;	

  //snapshot time and timestep
  float t;
  float dt;
  
};

typedef struct AthenaHeaderDef AthenaHeader;

/*! \fn AthenaHeader *ReadAthenaHeader(char *fname);
 *  \brief Function to read an Athena binary file header */
AthenaHeader *ReadAthenaHeader(FILE *fp);

/*! \fn WriteAthenaHeader(char *fname, AthenaHeader *h);
 *  \brief Function to write an Athena binary file header */
void WriteAthenaHeader(FILE *fp, AthenaHeader *h);

/*  \fn void ShowAthenaHeader(AthenaHeader *h)
 *  \brief Function to print an Athena Header to screen */
void ShowAthenaHeader(AthenaHeader *h);

#endif /* READ_ATHENA_HEADER_H */
