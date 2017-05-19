/*! \file read_athena_binary.h
 *  \brief Functions for reading athena binaries.
 *
 *  The individual athena binary fields (density, 
 *  velocity, etc.) are stored as grid_fft grids
 *  or fields, which eases use with Fourier 
 *  transform operations in k-space. This method
 *  is naturally parallelized using mpi.*/
#ifndef READ_ATHENA_TRACERS_H
#define READ_ATHENA_TRACERS_H

#define NO_PID
#define BAROTROPIC

typedef float Real;
typedef struct Tracer_s{
  Real d;                       /*!< density */
  Real M1;                      /*!< momentum density in 1-direction*/
  Real M2;                      /*!< momentum density in 2-direction*/
  Real M3;                      /*!< momentum density in 3-direction*/
#ifndef BAROTROPIC
  Real E;                       /*!< total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real B1c;                     /*!< cell centered magnetic fields in 1-dir*/
  Real B2c;                     /*!< cell centered magnetic fields in 2-dir*/
  Real B3c;                     /*!< cell centered magnetic fields in 3-dir*/
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];             /*!< passively advected scalars */
#endif
  Real x1;                      /*!< position in 1-direction*/
  Real x2;                      /*!< position in 2-direction*/
  Real x3;                      /*!< position in 3-direction*/
  long id;   /*!< tracer particle id */
#ifndef NO_PID
  int  pid;  /*!< orginal process id */
#endif /*NO_PID*/
}TracerS;

#endif /* READ_ATHENA_TRACERS_H */
