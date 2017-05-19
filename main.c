#include <stdio.h>
#include <random>
#include "kdtree2.hpp"
#include "write_shock_catalogues.hpp"

//typedef a boost 2-d array
typedef multi_array<float,2> array2dfloat;

int main(int argc, char **argv)
{
  int ir;
  int nr = 200;
  double r_min = 0.0;
  double r_max = 1.0;
  double *xi;
  double *r;
  double dr = (r_max-r_min)/((double) nr);
  double *DD;
  double *DR;
  double *RR;
  double *CR;
  double *CC;
  double rr;
  double rr_min = 1.0e-8;

  int imin = 0; // do we include autocorrelation of object w/ itself?
  int iminlim = 0;
  char fname[200];
  char fname_list[200];
  char fname_data[200];
  long nt = 0;
  vector<shock> s;
  vector<tracer> t;

  int ndim = 3;

  int n_R = 1000000;
  int n_D = 1000;

  kdtree2* tree_D;
  kdtree2* tree_R;
  kdtree2* tree_C;

  array2dfloat data_D;
  array2dfloat data_C;
  array2dfloat data_R;
  vector<float> qv(ndim);
  //result vectors
  kdtree2_result_vector res;

  xi  = (double *) malloc(nr*sizeof(double));
  r   = (double *) malloc(nr*sizeof(double));
  DD  = (double *) calloc(nr,sizeof(double));
  DR  = (double *) calloc(nr,sizeof(double));
  RR  = (double *) calloc(nr,sizeof(double));
  CR  = (double *) calloc(nr,sizeof(double));
  CC  = (double *) calloc(nr,sizeof(double));


  //initialize radial array
  for(int i=0;i<nr;i++)
  {
    r[i] = 0.5*dr + dr*i;
    //printf("i %d rmin %e r %e rmax %e\n",i,r[i]-0.5*dr,r[i],r[i]+0.5*dr);
  }

  std::default_random_engine rgen;
  std::uniform_real_distribution<double> dis(-1.0,2.0);
  std::uniform_real_distribution<double> disb(0.0,1.0);



  //print shock list filename
  sprintf(fname_list,"data/peak.blended.0753.list");
  sprintf(fname_data,"data/peak.blended.0753.dat");

  //read in shock list
  printf("Reading shock list %s.\n",fname_list);

  read_shock_list(fname_list, &s);
  printf("Number of shocks = %ld\n",s.size());
  for(int i=0;i<s.size();i++)
    nt += s[i].l;
  t.resize(nt);

  printf("Reading shock data %s.\n",fname_data);
  read_shock_data(fname_data, s, &t);

  printf("Generating data_D data....\n");
  fflush(stdout);

  n_D = 27*s.size();
  data_D.resize(extents[n_D][ndim]);

  //first, we have to duplicate the 
  //data by a factor of 27 -- this is
  //the stupid but easy to code way to do the
  //correlation function
  long j = 0;

  //0,1 first
  for(int i=0;i<s.size();i++)
  {
    data_D[j][0] = t[s[i].o].x[0];
    data_D[j][1] = t[s[i].o].x[1];
    data_D[j][2] = t[s[i].o].x[2];
    j++;
  }

  //replicate
  for(int ii=-1;ii<=1;ii++)
    for(int jj=-1;jj<=1;jj++)
      for(int kk=-1;kk<=1;kk++)
      {
        if(!((ii==0)&&(jj==0)&&(kk==0)))
          for(int i=0;i<s.size();i++)
          {
            data_D[j][0] = t[s[i].o].x[0] + ((double) ii);
            data_D[j][1] = t[s[i].o].x[1] + ((double) jj);
            data_D[j][2] = t[s[i].o].x[2] + ((double) kk);
            j++;
          }
      }
  printf("s->size %ld j %ld\n",s.size(),j);


  printf("Building data_D tree....\n");
  fflush(stdout);

  /*build the tree*/
  tree_D = new kdtree2(data_D);

  //generate random data
  printf("Generating R data....\n");
  fflush(stdout);

  n_R = n_D;

  data_R.resize(extents[n_R][ndim]);
  data_C.resize(extents[n_R][ndim]);

  //0,1 first
  j = 0;
  for(int i=0;i<s.size();i++)
  {
    for(int k=0;k<3;k++)
      data_R[j][k] = disb(rgen);
    for(int k=0;k<3;k++)
      data_C[j][k] = disb(rgen);    
    j++;
  }

  //replicate
  for(int ii=-1;ii<=1;ii++)
    for(int jj=-1;jj<=1;jj++)
      for(int kk=-1;kk<=1;kk++)
      {
        if(!((ii==0)&&(jj==0)&&(kk==0)))
          for(int i=0;i<s.size();i++)
          {
            data_R[j][0] = data_R[i][0] + ((double) ii);
            data_R[j][1] = data_R[i][1] + ((double) jj);
            data_R[j][2] = data_R[i][2] + ((double) kk);
            data_C[j][0] = data_C[i][0] + ((double) ii);
            data_C[j][1] = data_C[i][1] + ((double) jj);
            data_C[j][2] = data_C[i][2] + ((double) kk);
            j++;
          }
      }

  printf("Building data_R tree....\n");
  fflush(stdout);

  /*build the tree*/
  tree_R = new kdtree2(data_R);

  printf("Building data_C tree....\n");
  fflush(stdout);
  tree_C = new kdtree2(data_C);


  //correlation function

  //loop over objects
  for(long ip=0;ip<s.size();ip++)
  {
    printf("D ip %ld s.size() %ld\n",ip,s.size());
    for(int k=0;k<3;k++) 
      qv[k] = t[s[ip].o].x[k];
    //printf("r_max = %e\n",r_max);
    tree_D->r_nearest(qv, r_max, res);
    //printf("res.size() %ld\n",res.size());
    imin=iminlim;
    if(ip>0)
    {
      imin = 0;
    }
    for(long i=imin;i<res.size();i++)
    {
      rr = sqrt(res[i].dis);
      ir = (int) floor(rr/dr);
      //printf("i %ld res[i].dis %e ir %d\n",i,sqrt(res[i].dis),ir);
      if(rr>rr_min)
        DD[ir]+=1.0;
    }
  }

  //loop over objects
  for(long ip=0;ip<s.size();ip++)
  {
    printf("DR ip %ld s.size() %ld\n",ip,s.size());
    for(int k=0;k<3;k++) 
      qv[k] = t[s[ip].o].x[k];
    //printf("r_max = %e\n",r_max);
    tree_R->r_nearest(qv, r_max, res);
    //printf("res.size() %ld\n",res.size());
    imin=iminlim;
    if(ip>0)
    {
      imin = 0;
    }
    for(long i=imin;i<res.size();i++)
    {
      rr = sqrt(res[i].dis);
      ir = (int) floor(rr/dr);
      //printf("i %ld res[i].dis %e ir %d\n",i,sqrt(res[i].dis),ir);
      //never the same for cross correlation
      DR[ir]+=1.0;
    }
  }

  //loop over randoms
  for(long ip=0;ip<s.size();ip++)
  {
    printf("RR ip %ld s.size() %ld\n",ip,s.size());
    for(int k=0;k<3;k++) 
      qv[k] = data_R[ip][k];

    tree_R->r_nearest(qv, r_max, res);
    //printf("res.size() %ld\n",res.size());
    imin=iminlim;
    if(ip>0)
    {
      imin = 0;
    }
    for(long i=imin;i<res.size();i++)
    {
      rr = sqrt(res[i].dis);
      ir = (int) floor(rr/dr);
      //printf("i %ld res[i].dis %e ir %d\n",i,sqrt(res[i].dis),ir);
      if(rr>rr_min)
        RR[ir]+=1.0;
    }
  }
  //loop over randoms
  for(long ip=0;ip<s.size();ip++)
  {
    printf("CR ip %ld s.size() %ld\n",ip,s.size());
    for(int k=0;k<3;k++) 
      qv[k] = data_C[ip][k];

    tree_R->r_nearest(qv, r_max, res);
    //printf("res.size() %ld\n",res.size());
    imin=1;
    if(ip>0)
    {
      imin = iminlim;
    }
    for(long i=imin;i<res.size();i++)
    {
      rr = sqrt(res[i].dis);
      ir = (int) floor(rr/dr);
      //printf("i %ld res[i].dis %e ir %d\n",i,sqrt(res[i].dis),ir);
      CR[ir]+=1.0;
    }
  }

  //loop over randoms
  for(long ip=0;ip<s.size();ip++)
  {
    printf("CC ip %ld s.size() %ld\n",ip,s.size());
    for(int k=0;k<3;k++) 
      qv[k] = data_C[ip][k];

    tree_C->r_nearest(qv, r_max, res);
    //printf("res.size() %ld\n",res.size());
    imin=iminlim;
    if(ip>0)
    {
      imin = 0;
    }
    for(long i=imin;i<res.size();i++)
    {
      rr = sqrt(res[i].dis);
      ir = (int) floor(rr/dr);
      //printf("i %ld res[i].dis %e ir %d\n",i,sqrt(res[i].dis),ir);
      if(rr>rr_min)
        CC[ir]+=1.0;
    }
  }
  //print DD, DR, RR
  for(int i=0;i<nr;i++)
    printf("i %d r %e DD %e DR %e RR %e RC %e CC %e\n",i,r[i],DD[i],DR[i],RR[i],CR[i],CC[i]);

  //save data
  sprintf(fname,"correlation.txt");
  FILE *fp;
  fp = fopen(fname,"w");
  for(int i=0;i<nr;i++)
    fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\n",r[i],DD[i],DR[i],RR[i],CR[i],CC[i]);
  fclose(fp);

  //free memory
  free(xi);
  free(r);
  free(DD);
  free(DR);
  free(RR);
  free(CR);
  free(CC);



  return 0;
}