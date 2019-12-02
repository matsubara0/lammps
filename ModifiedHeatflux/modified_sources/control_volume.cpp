/* -*- c++ -*- ----------------------------------------------------------
LAST_MODIFIED="2018/12/11 11:09:45" 
  
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

#include <math.h>
#include <string.h>
#include "control_volume.h"
#include "comm.h"
#include "domain.h"
#include "region.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
//using namespace MathConst;

/* ---------------------------------------------------------------------- */

ControlVolume::ControlVolume(LAMMPS *lmp, int max_cv_regions_in) : Pointers(lmp)
{
  allocated = 0;
  max_cv_regions = max_cv_regions_in;
  ncvs = 0;
  iregion = NULL;
  idregion = NULL;
  cvzm = NULL;
  cvzp = NULL;
  allocate();
}

ControlVolume::~ControlVolume()
{
  if (allocated) {
    memory->destroy(iregion);
    memory->destroy(idregion);
    memory->destroy(cvzm);
    memory->destroy(cvzp);
  }
}

void ControlVolume::allocate()
{
  allocated = 1;
  int l_str = 5;
  memory->create(iregion ,max_cv_regions,"pair:iregion");
  memory->create(idregion,max_cv_regions,l_str,"pair:idregion");
  memory->create(cvzm,max_cv_regions,"pair:cvzm");
  memory->create(cvzp,max_cv_regions,"pair:cvzp");
}

int ControlVolume::register_cv(char* region_name)
{
  int rid = domain->find_region(region_name);

  for(int i=0;i<ncvs; i++){
    if(rid == iregion[i]) return i; //skip if already defined
  }

  int icv = ncvs;
  ncvs++;
  if(ncvs > max_cv_regions)  error->all(FLERR,"The maxmum number of CVs is exceeded");

  iregion[icv] = rid;
  if(rid <=-1) error->all(FLERR,"Region ID for CV does not exist");
  strcpy(idregion[icv],region_name);

  Region *region = NULL;
  region = domain->regions[rid];
  cvzm[icv] = region->extent_zlo;
  cvzp[icv] = region->extent_zhi;

  FILE *out = screen;
  if(comm->me == 0){
    for(int d=0;d<2;d++){//2 for screen and logfile
      if(out){
	fprintf(out,"CV%d defined: region no=%d, id=%s, range %g <= z <=%g\n"
		,icv, iregion[icv],idregion[icv], cvzm[icv], cvzp[icv]);
      }
      out = logfile;
    }
  }

  return icv;
  
}

double ControlVolume::get_rfraction(int i, double zi, double zj)
{
  double frac=0.0; //fraction of rij must be 0<= frac <=1

  if(iregion[i]<0) return 1.0; //CV is whole space

  if(zi < cvzm[i]){
    if(zj < cvzm[i]){                // ij| |
      frac = 0.0;
    }else if( zj <= cvzp[i]){        // i |j|
      frac = (zj - cvzm[i])/(zj-zi);
      }else{                           // i | |j
      frac = (cvzp[i]-cvzm[i])/(zj-zi);
    }
  }else if(zi <= cvzp[i]){
    if(zj < cvzm[i]){                // j|i |
      frac = (zi-cvzm[i])/(zi-zj);
    }else if(zj <= cvzp[i]){         //  |ij|
      frac = 1.0;
    }else{                           //  |i |j
      frac = (cvzp[i]-zi)/(zj-zi);
    }
  }else{
    if(zj < cvzm[i]){                // j| |i
      frac = (cvzp[i]-cvzm[i])/(zi-zj);
    }else if(zj <= cvzp[i]){         //  |j|i
      frac = (cvzp[i]-zj)/(zi-zj);
    }else{                           //  | |ij
      frac = 0.0;
    }
  }
  if (frac < 0.0) error->all(FLERR,"frac is negative.");
  return frac;
}

//double ControlVolume::memory_usage()
//{
//  double bytes = comm->nthreads*ncvs *sizeof(int); //iregion
//  bytes += comm->nthreads*ncvs*5 * sizeof(double);  //idregion
//  bytes += comm->nthreads*ncvs*2 * sizeof(double); //cvzm,cvzp
//  return bytes;
//}
