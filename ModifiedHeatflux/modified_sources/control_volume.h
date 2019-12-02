/* -*- c++ -*- ----------------------------------------------------------
LAST_MODIFIED="2018/12/10 19:12:42" 
  
   H.M. 

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef LMP_CONTROL_VOLUME_H
#define LMP_CONTROL_VOLUME_H

#include <math.h>
#include "pointers.h"

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */
class ControlVolume : protected Pointers {

 public:
  int allocated;   //if array for CV regions are allocated
  int ncvs;   // number of CV
  int *iregion;    //region ids for CVs
  char **idregion; //name of CV regions
  double *cvzm;     //left z-coordinates CVs
  double *cvzp;     //right z-coordinates CVs

  ControlVolume(class LAMMPS *, int max_cv_regions_in);
  virtual ~ControlVolume();
  virtual void allocate();
  //  double memory_usage();
  int register_cv(char* region_name);
  double get_rfraction(int i, double zi, double zj);
  
  inline int size_vatom_cols(){return 9*ncvs;}
  
private:
  int max_cv_regions;

};

}

#endif
