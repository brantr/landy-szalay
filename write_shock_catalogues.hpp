#ifndef  SHOCK_CATALOGUES_H
#define  SHOCK_CATALOGUES_H

#include"shock_data_types.hpp"


void read_tracer_maxima(char fdir_out[], int isnap, int isub, float *x_min, float *x_max);

long load_tracers(char fdir[], char filebase[], char fsuffix[], char fdir_out[], int isnap, int isub, vector<tracer> *t, float dthresh);

void read_shock_list_isub(char fdir_out[], int isnap, int isub, vector<shock> *s);
void read_shock_data_isub(char fdir_out[], int isnap, int isub, vector<shock> s, vector<tracer> *t);
void read_shock_list_isnap(char fdir_out[], int isnap, vector<shock> *s);
void read_shock_data_isnap(char fdir_out[], int isnap, vector<shock> s, vector<tracer> *t);
void read_shock_list(char fname[], vector<shock> *s);
void read_shock_data(char fname[], vector<shock> s, vector<tracer> *t);

void write_shock_list_isub(char fdir_out[], int isnap, int isub, vector<shock> s);
void write_shock_data_isub(char fdir_out[], int isnap, int isub, vector<shock> s, vector<tracer> t);
void write_shock_list_isnap(char fdir_out[], int isnap, vector<shock> s);
void write_shock_data_isnap(char fdir_out[], int isnap, vector<shock> s, vector<tracer> t);
void write_shock_list(char fname[], vector<shock> s);
void write_shock_data(char fname[], vector<shock> s, vector<tracer> t);

void write_null_shock_list_isnap(char fdir_out[], int isnap);
void write_null_shock_data_isnap(char fdir_out[], int isnap);
void write_null_shock_list_isub(char fdir_out[], int isnap, int isub);
void write_null_shock_data_isub(char fdir_out[], int isnap, int isub);
void write_null_shock_list_nfiles(char fdir_out[], int isnap, int nfiles);
void write_null_shock_data_nfiles(char fdir_out[], int isnap, int nfiles);
void write_null_shock_list(char fname[]);
void write_null_shock_data(char fname[]);

void write_merge_list(char fdir_out[], int n, vector<merger> m);
int read_merge_list(char fdir_out[], vector<merger> *m);

void load_ordered_shock_data(char fdir_out[], int isnap, vector<shock> s, vector<shock_sort> ss, vector<tracer> *t);


#endif /*SHOCK_CATALOGUES_H*/
