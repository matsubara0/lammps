Jan. 12, 2021 Hiroki Matsubara

This folder includes the modified source codes of lammps to compute heat flux vector that satisfies energy conservation. For details, please refer to : H. Matsubara, G. Kikugawa, T. Bessho, T. Ohara, Diamond & Related Materials (under review).

** Please note that the current version is not designed for general use.
(For general use, in future, we are planning to provide a modified version of pair/airebo that is compatible with the centroid approximation centroid/stress/atom. See https://lammps.sandia.gov/doc/compute_stress_atom.html)

1) Contents
- The directory modified_sources includes the modified codes.
- The directory original_sources includes the original lammps codes (of lammps-7Aug19) that are necessary to recover the original lammps.
- The directory sample includes the input file that reproduce results of the reference paper above for diamond crystal.

2) Disclaimer of Warranty
Use of this code is at your sole risk. Anyone can use these codes, but you must be careful because the current version of the modified codes are not designed for general use. Source codes are provided "as is" with no warranties or guarantees.

3) How to Use

3-1)Our code assumes to be used with the lammps version lammps-7Aug19, and supports pair/airebo, pair/tersoff, and pair/lj/cut only. Use of any other potential type (including angle, dihedral, and pair/hybrid etc.) can cause a significant error.

3-2) In the input file, use compute stress/atom/local insted of compute stress/atom as below.

 region cv block INF INF INF INF ${z1} ${z2}
 group cvatom dynamic all region cv every 1
 compute mystress all stress/atom/local cv NULL pair

Here, we assume that the average heat flux of a control volume 'cv' is mesasrued. Only one such control volume can be defined in a single lammps input file.
The group for stress/atom/local should be "all" because atoms outside cv can contribute to the heat flux if their interactions cross cv.

3-3) The heat flux can be measured in all directions, but the control volume can be finite only in the z direction.

3-4)Please do not use 'compute stress/atom' if 'compute stress/atom/local' is being used.

3-5)To install the codes, just copy the source (*.cpp) and header (*.h) files in the lammps source directory and replace tersoff.* and airebo.* in MANYBODY package as:

 cp modified_source/{*.cpp,*.h} lammps-7Aug19/src
 cp modified_source/{tersoff.*,airebo.*} lammps-7Aug19/src/MANYBODY
 
Then, make as usual. Note that this will overwirte some of the original lammps codes. The following commands recovers the original lammps codes.
 
 cp original_source/{*.cpp,*.h} lammps-7Aug19/src
 cp origianl_source/{tersoff.*,airebo.*} lammps-7Aug19/src/MANYBODY
 rm lammps-7Aug19/src/{compute_stress_atom_local.*,control_volume.*}
