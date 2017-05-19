#include <stdio.h>
#include "write_shock_catalogues.hpp"
#include "read_athena_header.hpp"
#include "read_athena_tracers.h"

#ifdef SUBSAMPLE
#include "rng.h"
#endif /*SUBSAMPLE*/

using namespace std;

void write_merge_list(char fdir_out[], int n, vector<merger> m)
{
  //write the merge list to file
  FILE *fp;
  char fmerge[200];
  sprintf(fmerge,"%s/merge.dat",fdir_out);
  if(!(fp = fopen(fmerge,"w")))
  {
    printf("Error opening %s\n",fmerge);
    fflush(stdout);
    exit(-1); 
  }
   
  fwrite(&n,1,sizeof(int),fp);
  fwrite(&m[0],n,sizeof(merger),fp);
  fclose(fp);
}
int read_merge_list(char fdir_out[], vector<merger> *m)
{
  //read the merge from file
  FILE *fp;
  char fmerge[200];
  int n;
  sprintf(fmerge,"%s/merge.dat",fdir_out);
  merger min;
  if(!(fp = fopen(fmerge,"r")))
  {
    printf("Error opening %s\n",fmerge);
    fflush(stdout);
    exit(-1); 
  }
   
  fread(&n,1,sizeof(int),fp);
  m->resize(n);
  for(int i=0;i<n;i++)
  {
    fread(&min,1,sizeof(merger),fp);
    (*m)[i] = min;
  }
  fclose(fp);

  return n;
}


void read_tracer_maxima(char fdir_out[], int isnap, int isub, float *x_min, float *x_max)
{
  FILE *fp;
  char filename[200];

  //write min and max information to output directory
  sprintf(filename,"%s/range.%04d.%04d.txt",fdir_out,isnap,isub);
  if(!(fp = fopen(filename,"r")))
  {
    printf("Error opening %s.\n",filename);
    fflush(stdout);
  }
  fscanf(fp,"%f\t%f\t%f\n",&x_min[0],&x_min[1],&x_min[2]);
  fscanf(fp,"%f\t%f\t%f\n",&x_max[0],&x_max[1],&x_max[2]);
  fclose(fp);
}

long load_tracers(char fdir[], char filebase[], char fsuffix[], char fdir_out[], int isnap, int isub, vector<tracer> *t, float dthresh)
{
  FILE *fp;
  char filename[200];   /*name of file containing the tracers*/
  long n_tracers;       /*number of tracers in the file*/
  float *d;             /*tracer densities*/
  float *x;             /*tracer x positions*/
  float *y;             /*tracer y positions*/
  float *z;             /*tracer z positions*/
  float *vx;            /*tracer x velocities*/
  float *vy;            /*tracer y velocities*/
  float *vz;            /*tracer z velocities*/
  long  *l;             /*tracer ids*/
  long ntd = 0;         /*number of tracers above the density threshold*/
  AthenaHeader *h;

  //corners of the bounding box
  float t_min[3] = {1e9,1e9,1e9};
  float t_max[3] = {-1e9,-1e9,-1e9};


  //buffer for storing tracers into the tracer vector *t
  tracer tin;
  
  if(isub==0)
  {
    /*create a new filename*/
    sprintf(filename,"%s/%s.%04d.%s",fdir,filebase,isnap,fsuffix);
  }else{
    /*create a new filename*/
    sprintf(filename,"%s/%s-id%d.%04d.%s",fdir,filebase,isub,isnap,fsuffix);
  }

  /*open tracer file*/
  if(!(fp = fopen(filename,"r")))
  {
    printf("Error opening %s in load tracers (fdir=%s, filebase=%s, fdout=%s.\n",filename,fdir,filebase,fdir_out);
    fflush(stdout); 
  }else{
    printf("Opening %s.\n",filename);
    fflush(stdout);
  }

  /* Read header */
  h = ReadAthenaHeader(fp);

  ShowAthenaHeader(h);


  /* read the number of tracer in this file */ 
  fread(&n_tracers,1,sizeof(long),fp);

  printf("n_tracers = %ld\n",n_tracers);
  fflush(stdout);

  /* Allocate buffer */
  if(!(d = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property d (n_tracers = %ld).\n",n_tracers);
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(x = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property buf.\n");
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(y = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property y.\n");
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(z = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property z.\n");
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(vx = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property z.\n");
    fflush(stdout);
    exit(-1);
  }

    /* Allocate buffer */
    if(!(vy = (float *)malloc(n_tracers*sizeof(float))))
    {
      printf("Error allocating tracer property z.\n");
      fflush(stdout);
      exit(-1);
    }


    /* Allocate buffer */
    if(!(vz = (float *)malloc(n_tracers*sizeof(float))))
    {
      printf("Error allocating tracer property z.\n");
      fflush(stdout);
      exit(-1);
    }


    /* read density */
    fread(d,n_tracers,sizeof(float),fp);

    /* read M1 */
    fread(x,n_tracers,sizeof(float),fp);
    for(long tt=0;tt<n_tracers;tt++)
      vx[tt] = x[tt]/d[tt];
  
    /* read M2 */
    fread(x,n_tracers,sizeof(float),fp);
    for(long tt=0;tt<n_tracers;tt++)
      vy[tt] = x[tt]/d[tt];

    /* read M3 */
    fread(x,n_tracers,sizeof(float),fp);
    for(long tt=0;tt<n_tracers;tt++)
      vz[tt] = x[tt]/d[tt];

#ifndef BAROTROPIC
    /* read E */
    fread(x,n_tracers,sizeof(float),fp);
#endif /* BAROTROPIC */

#ifdef MHD
    /* read B1c */
    fread(x,n_tracers,sizeof(float),fp);

    /* read B2c */
    fread(x,n_tracers,sizeof(float),fp);

    /* read B3c */
    fread(x,n_tracers,sizeof(float),fp);
#endif /*MHD*/

#if (NSCALARS > 0)
    for(k=0;k<NSCALARS;k++)
      fread(x,n_tracers,sizeof(float),fp);
#endif

    /* read x1 */
    fread(x,n_tracers,sizeof(float),fp);
  
    /* read x2 */
    fread(y,n_tracers,sizeof(float),fp);

    /* read x3 */
    fread(z,n_tracers,sizeof(float),fp);  

    /* Allocate buffer */
    if(!(l = (long *)malloc(n_tracers*sizeof(long))))
    {
      printf("Error allocating tracer property buf.\n");
      fflush(stdout);
    }

    /* read id */
    fread(l,n_tracers,sizeof(long),fp);

    /*close tracer file*/
    fclose(fp);


    /*keep only particles with densities above or = threshold*/
    ntd = 0;
    for(long tt=0;tt<n_tracers;tt++)
    {

      if( x[tt] < t_min[0])
        t_min[0] = x[tt];
      if( y[tt] < t_min[1])
        t_min[1] = y[tt];
      if( z[tt] < t_min[2])
        t_min[2] = z[tt];

      if( x[tt] > t_max[0])
        t_max[0] = x[tt];
      if( y[tt] > t_max[1])
        t_max[1] = y[tt];
      if( z[tt] > t_max[2])
        t_max[2] = z[tt];

      //if particle is above the threshold, keep it
      if(d[tt]>=dthresh)
      {  
        tin.id = l[tt];
        tin.d = d[tt];
        tin.x[0] = x[tt];
        tin.x[1] = y[tt];
        tin.x[2] = z[tt];
        tin.v[0] = vx[tt];
        tin.v[1] = vy[tt];
        tin.v[2] = vz[tt];

#ifdef SUBSAMPLE
      //if(rng_uniform(0,1)<0.125)
        //if(rng_uniform(0,1)<20)
        if(1)
        //if(rng_uniform(0,1)<1./64.)
        {
#endif /*SUBSAMPLE*/

        //add to tracer list
        (*t).push_back(tin);

        //remember that we've kept a particle
        ntd++;
#ifdef SUBSAMPLE
        }//rng_uniform
#endif /*SUBSAMPLE*/

      }
    }

    /*free buffer memory*/
    free(d);
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(l);


  /*free header*/
  free(h);


  //write min and max information to output directory
  sprintf(filename,"%s/range.%04d.%04d.txt",fdir_out,isnap,isub);
  if(!(fp = fopen(filename,"w")))
  {
    printf("Error opening %s.\n",filename);
    fflush(stdout);
  }
  fprintf(fp,"%e\t%e\t%e\n",t_min[0],t_min[1],t_min[2]);
  fprintf(fp,"%e\t%e\t%e\n",t_max[0],t_max[1],t_max[2]);
  fclose(fp);

  //return number of tracers > density
  return ntd;

}

void read_shock_data_isub(char fdir_out[], int isnap, int isub, vector<shock> s, vector<tracer> *t)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.%04d.dat",fdir_out,isnap,isub);
  read_shock_data(fname,s,t);
}

void read_shock_data(char fname[], vector<shock> s, vector<tracer> *t)
{
  int  n;
  long io;
  FILE *fp_dat;


  float *dpeak;
  float *xpeak;
  float *ypeak;
  float *zpeak;
  float *vxpeak;
  float *vypeak;
  float *vzpeak;
  long  *lpeak;
  
  //write snap peak details

  if(!(fp_dat = fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    fflush(stdout);
  }
  fread(&n,sizeof(int),1,fp_dat);


  io = 0;
  for(int i=0;i<n;i++)
  {
    dpeak  = (float *) malloc(s[i].l * sizeof(float));
    xpeak  = (float *) malloc(s[i].l * sizeof(float));
    ypeak  = (float *) malloc(s[i].l * sizeof(float));
    zpeak  = (float *) malloc(s[i].l * sizeof(float));
    vxpeak = (float *) malloc(s[i].l * sizeof(float));
    vypeak = (float *) malloc(s[i].l * sizeof(float));
    vzpeak = (float *) malloc(s[i].l * sizeof(float));
    lpeak  = (long  *) malloc(s[i].l * sizeof(long));


    fread(dpeak,sizeof(float),s[i].l,fp_dat);
    fread(xpeak,sizeof(float),s[i].l,fp_dat);
    fread(ypeak,sizeof(float),s[i].l,fp_dat);
    fread(zpeak,sizeof(float),s[i].l,fp_dat);
    fread(vxpeak,sizeof(float),s[i].l,fp_dat);
    fread(vypeak,sizeof(float),s[i].l,fp_dat);
    fread(vzpeak,sizeof(float),s[i].l,fp_dat);
    fread(lpeak,sizeof(long),s[i].l,fp_dat);

    for(long k=0;k<s[i].l;k++)
    {
      (*t)[io].d    = dpeak[k];
      (*t)[io].x[0] = xpeak[k];
      (*t)[io].x[1] = ypeak[k];
      (*t)[io].x[2] = zpeak[k];
      (*t)[io].v[0] = vxpeak[k];
      (*t)[io].v[1] = vypeak[k];
      (*t)[io].v[2] = vzpeak[k];
      (*t)[io].id   = lpeak[k];

      io++;
    }

    free(dpeak);
    free(xpeak);
    free(ypeak);
    free(zpeak);
    free(vxpeak);
    free(vypeak);
    free(vzpeak);
    free(lpeak);
  }

  //close peak details
  fclose(fp_dat);
}

void read_shock_list_isub(char fdir_out[], int isnap, int isub, vector<shock> *s)
{

  char fname[200];
  sprintf(fname,"%s/peak.%04d.%04d.list",fdir_out,isnap,isub);  
  read_shock_list(fname,s);

}

void read_shock_list(char fname[], vector<shock> *s)
{

  FILE  *fp_list;
  int    n;
  long  *l;
  long  *o;
  float *d;
  long  *id;
  float *b_min;
  float *b_max;

  //write the peak list to file
  if(!(fp_list = fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    fflush(stdout);
  }

  //read number of shocks
  fread(&n,1,sizeof(int),fp_list);

  //resize s vector
  s->resize(n);

  l	= (long  *) malloc(n*sizeof(long));
  o     = (long  *) malloc(n*sizeof(long));
  d     = (float *) malloc(n*sizeof(float));
  id    = (long  *) malloc(n*sizeof(long));
  b_min = (float *) malloc(3*n*sizeof(float));
  b_max = (float *) malloc(3*n*sizeof(float));


  //read peak catalog for the B snapshot 
  fread(l,n,sizeof(long),fp_list);
  fread(o,n,sizeof(long),fp_list);
  fread(d,n,sizeof(float),fp_list);
  fread(id,n,sizeof(long),fp_list);
  fread(b_min,3*n,sizeof(float),fp_list);
  fread(b_max,3*n,sizeof(float),fp_list);


  for(int k=0;k<n;k++)
  {
    (*s)[k].l  = l[k];
    (*s)[k].o  = o[k];
    (*s)[k].d  = d[k];
    (*s)[k].id = id[k];

    for(int j=0;j<3;j++)
    {
      (*s)[k].min[j] = b_min[3*k+j];
      (*s)[k].max[j] = b_max[3*k+j];
    }
  }

  //close the file
  fclose(fp_list);

  //free memory
  free(l);
  free(o);
  free(d);
  free(id);
  free(b_min);
  free(b_max);
}

void write_null_shock_data_isnap(char fdir_out[], int isnap)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.dat",fdir_out,isnap);
  write_null_shock_data(fname);
}
void write_null_shock_data_isub(char fdir_out[], int isnap, int isub)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.%04d.dat",fdir_out,isnap,isub);
  write_null_shock_data(fname);
}
void write_null_shock_data_nfiles(char fdir_out[], int isnap, int nfiles)
{
  int isub;
  for(isub=0;isub<nfiles;isub++)
  {
    //write snap peak details
    write_null_shock_data_isub(fdir_out, isnap, isub);
  }

}

void write_null_shock_data(char fname[])
{
  int  n;
  FILE *fp_dat;

  if(!(fp_dat = fopen(fname,"w")))
  {
    printf("Error opening %s.\n",fname);
    fflush(stdout);
  }
  n = 0;
  fwrite(&n,sizeof(int),1,fp_dat);
  fclose(fp_dat);
}

void write_shock_data_isub(char fdir_out[], int isnap, int isub, vector<shock> s, vector<tracer> t)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.%04d.dat",fdir_out,isnap,isub);
  write_shock_data(fname,s,t);
}

void write_shock_data_isnap(char fdir_out[], int isnap, vector<shock> s, vector<tracer> t)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.dat",fdir_out,isnap);
  write_shock_data(fname,s,t);
}

void write_shock_data(char fname[], vector<shock> s, vector<tracer> t)
{
  int  n;
  long io;
  FILE *fp_dat;

  float *dpeak;
  float *xpeak;
  float *ypeak;
  float *zpeak;
  float *vxpeak;
  float *vypeak;
  float *vzpeak;
  long  *lpeak;
  
  //write snap peak details

  if(!(fp_dat = fopen(fname,"w")))
  {
    printf("Error opening %s.\n",fname);
    fflush(stdout);
  }
  io = 0;
  n = s.size();

  //printf("s.size() %ld\n",s.size());
  //fflush(stdout);

  fwrite(&n,sizeof(int),1,fp_dat);
  if(n==0)
  {
    fclose(fp_dat);
    return;
  }

  for(int i=0;i<n;i++)
  {
    dpeak  = (float *) malloc(s[i].l * sizeof(float));
    xpeak  = (float *) malloc(s[i].l * sizeof(float));
    ypeak  = (float *) malloc(s[i].l * sizeof(float));
    zpeak  = (float *) malloc(s[i].l * sizeof(float));
    vxpeak = (float *) malloc(s[i].l * sizeof(float));
    vypeak = (float *) malloc(s[i].l * sizeof(float));
    vzpeak = (float *) malloc(s[i].l * sizeof(float));
    lpeak  = (long  *) malloc(s[i].l * sizeof(long));

    for(long k=0;k<s[i].l;k++)
    {
      dpeak[k]  = t[io].d;
      xpeak[k]  = t[io].x[0];
      ypeak[k]  = t[io].x[1];
      zpeak[k]  = t[io].x[2];
      vxpeak[k] = t[io].v[0];
      vypeak[k] = t[io].v[1];
      vzpeak[k] = t[io].v[2];
      lpeak[k]  = t[io].id;

		
      io++;
    }

    fwrite(dpeak,sizeof(float),s[i].l,fp_dat);
    fwrite(xpeak,sizeof(float),s[i].l,fp_dat);
    fwrite(ypeak,sizeof(float),s[i].l,fp_dat);
    fwrite(zpeak,sizeof(float),s[i].l,fp_dat);
    fwrite(vxpeak,sizeof(float),s[i].l,fp_dat);
    fwrite(vypeak,sizeof(float),s[i].l,fp_dat);
    fwrite(vzpeak,sizeof(float),s[i].l,fp_dat);
    fwrite(lpeak,sizeof(long),s[i].l,fp_dat);

    free(dpeak);
    free(xpeak);
    free(ypeak);
    free(zpeak);
    free(vxpeak);
    free(vypeak);
    free(vzpeak);
    free(lpeak);
  }

  //close peak details
  fclose(fp_dat);
}

void write_null_shock_list_nfiles(char fdir_out[], int isnap, int nfiles)
{
  int isub;

  for(isub=0;isub<nfiles;isub++)
    write_null_shock_list_isub(fdir_out, isnap, isub);
}

void write_null_shock_list_isub(char fdir_out[], int isnap, int isub)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.%04d.list",fdir_out,isnap,isub);
  write_null_shock_list(fname);
}

void write_null_shock_list_isnap(char fdir_out[], int isnap)
{
  char fname[200];
  sprintf(fname,"%s/peak.%04d.list",fdir_out,isnap);
  write_null_shock_list(fname);
}

void write_null_shock_list(char fname[])
{
  FILE  *fp_list;
  int    n;
  long  *l;
  long  *o;
  float *d;
  long  *id;
  float *b_min;
  float *b_max;

  //write the peak list to file

  if(!(fp_list = fopen(fname,"w")))
  {
    printf("Error opening %s.\n",fname);
    fflush(stdout);
  }

  n = 0;

#if DEBUG >= 3
  printf("WRITE_NULL_LIST %s\n",fname);
  fflush(stdout);
#endif //DEBUG

  fwrite(&n,1,sizeof(int),fp_list);
  fclose(fp_list);
}


void write_shock_list_isub(char fdir_out[], int isnap, int isub, vector<shock> s)
{
  //file name string
  char fname[200];

  //make name for file
  sprintf(fname,"%s/peak.%04d.%04d.list",fdir_out,isnap,isub);

  //write the shock list
  write_shock_list(fname,s);
}

void write_shock_list_isnap(char fdir_out[], int isnap, vector<shock> s)
{
  //file name string
  char fname[200];

  //make name for file
  sprintf(fname,"%s/peak.%04d.list",fdir_out,isnap);

  //write the shock list
  write_shock_list(fname,s);
}

void write_shock_list(char fname[], vector<shock> s)
{
  FILE  *fp_list;
  int   n;
  long  *l;
  long  *o;
  float *d;
  long  *id;
  float *b_min;
  float *b_max;

  //open the file
  if(!(fp_list = fopen(fname,"w")))
  {
    printf("Error opening %s.\n",fname);
    fflush(stdout);
  }

  n	= (int) s.size();

  fwrite(&n,1,sizeof(int),fp_list);
  if(n==0)
  {
    fclose(fp_list);  
    return;
  }


  l	    = (long  *) malloc(n*sizeof(long));
  o     = (long  *) malloc(n*sizeof(long));
  d     = (float *) malloc(n*sizeof(float));
  id    = (long  *) malloc(n*sizeof(long));
  b_min = (float *) malloc(3*n*sizeof(float));
  b_max = (float *) malloc(3*n*sizeof(float));

  for(int k=0;k<n;k++)
  {
    l[k]  = s[k].l;
    o[k]  = s[k].o;
    d[k]  = s[k].d;
    id[k] = s[k].id;

#if DEBUG >= 3
    printf("WRITE_SHOCK_LIST k %3d l %10ld o %10ld d %9.8e id %10ld %s\n",k,l[k],o[k],d[k],s[k].id,fname);
    fflush(stdout);
#endif

    for(int j=0;j<3;j++)
    {
      b_min[3*k+j] = s[k].min[j];
      b_max[3*k+j] = s[k].max[j];
    }
  }

  //read peak catalog for the B snapshot 
  fwrite(l,n,sizeof(long),fp_list);
  fwrite(o,n,sizeof(long),fp_list);
  fwrite(d,n,sizeof(float),fp_list);
  fwrite(id,n,sizeof(long),fp_list);
  fwrite(b_min,3*n,sizeof(float),fp_list);
  fwrite(b_max,3*n,sizeof(float),fp_list);

  //close the file
  fclose(fp_list);

  //free memory
  free(l);
  free(o);
  free(d);
  free(id);
  free(b_min);
  free(b_max);
}




void load_ordered_shock_data(char fdir_out[], int isnap, vector<shock> s, vector<shock_sort> ss, vector<tracer> *t)
{
  int  n;
  long io = 0;
  FILE *fp_dat;
  char fname[200];

  float *dpeak;
  float *xpeak;
  float *ypeak;
  float *zpeak;
  float *vxpeak;
  float *vypeak;
  float *vzpeak;
  long  *lpeak;


  int isub;
  long ishock;

  long offset;


  
  for(ishock=0;ishock<s.size();ishock++)
  {
    isub = ss[ishock].isub;

  
    //write snap peak details
    sprintf(fname,"%s/peak.%04d.%04d.dat",fdir_out,isnap,isub);
    if(!(fp_dat = fopen(fname,"r")))
    {
      printf("Error opening %s in load_ordered_shock_data.\n",fname);
      fflush(stdout);
    }


    //move to an offset
    
    offset = sizeof(int) + s[ishock].o * (7*sizeof(float) + sizeof(long));
    fseek(fp_dat, offset, SEEK_SET);
  
    dpeak  = (float *) malloc(s[ishock].l * sizeof(float));
    xpeak  = (float *) malloc(s[ishock].l * sizeof(float));
    ypeak  = (float *) malloc(s[ishock].l * sizeof(float));
    zpeak  = (float *) malloc(s[ishock].l * sizeof(float));
    vxpeak = (float *) malloc(s[ishock].l * sizeof(float));
    vypeak = (float *) malloc(s[ishock].l * sizeof(float));
    vzpeak = (float *) malloc(s[ishock].l * sizeof(float));
    lpeak  = (long  *) malloc(s[ishock].l * sizeof(long));


    fread(dpeak, sizeof(float),s[ishock].l,fp_dat);
    fread(xpeak, sizeof(float),s[ishock].l,fp_dat);
    fread(ypeak, sizeof(float),s[ishock].l,fp_dat);
    fread(zpeak, sizeof(float),s[ishock].l,fp_dat);
    fread(vxpeak,sizeof(float),s[ishock].l,fp_dat);
    fread(vypeak,sizeof(float),s[ishock].l,fp_dat);
    fread(vzpeak,sizeof(float),s[ishock].l,fp_dat);
    fread(lpeak, sizeof(long), s[ishock].l,fp_dat);

    for(long k=0;k<s[ishock].l;k++)
    {
      (*t)[io].d    = dpeak[k];
      (*t)[io].x[0] = xpeak[k];
      (*t)[io].x[1] = ypeak[k];
      (*t)[io].x[2] = zpeak[k];
      (*t)[io].v[0] = vxpeak[k];
      (*t)[io].v[1] = vypeak[k];
      (*t)[io].v[2] = vzpeak[k];
      (*t)[io].id   = lpeak[k];

      io++;
    }

    free(dpeak);
    free(xpeak);
    free(ypeak);
    free(zpeak);
    free(vxpeak);
    free(vypeak);
    free(vzpeak);
    free(lpeak);

    //close peak details
    fclose(fp_dat);
  }
}

