#ifndef SHOCK_DATA_TYPES_H
#define SHOCK_DATA_TYPES_H
#include "kdtree2.hpp"

struct shock_sort 
{
  int isub;
  float d;
};

struct pos
{
  float x[3];
};

struct tracer
{
  long id;
  float d;
  float x[3];
  float v[3];
  long peak_index;
  //long peak_id;
  //long neighbor_index;
};

struct merger
{
  int snap_A;
  int i_A;
  float d_A;
  long  id_A;
  int snap_B;
  int i_B;
  float d_B;
  long id_B;
};

struct shock
{
  long l;   //length
  long o;   //offset
  float d;  //density
  long id;  //peak id, inherited from densest particle
  float min[3];	//bounding box min
  float max[3];	//bounding box max
};

struct box
{
  float min[3];
  float max[3];
};




// define, for convenience a 2d array of floats. 
//
typedef multi_array<float,2> array2dfloat;

//std::set<long> idc;
extern std::vector<long> idc;
extern std::vector<tracer> tv;
extern std::vector<tracer> td;
extern tracer tin;

extern kdtree2*	tree;
extern array2dfloat    data;

#endif /*SHOCK_DATA_TYPES_H*/
