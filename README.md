MFMM_tstar
Introduction

MFMM_tstar is a software package for computing the seismic attenuation operator t* using a modified fast marching method (MFMM).

Traditional methods for calculating t* require explicit ray tracing followed by integration along ray paths. In contrast, MFMM_tstar accumulates the attenuation operator directly during the evolution of the wavefront within the fast marching framework. This approach avoids explicit ray tracing and ensures consistent accumulation of attenuation along minimum‐traveltime paths.

This software is designed for studies of seismic attenuation and wave propagation in heterogeneous absorbing media.




MFMM_tstar

Introduction of 'MFMM_tstar'

tstar_FMM is a software for calculating attenuation operator t* directly by the fast marching method without ray tracing.
*** Please cite the software as:*** Dongdong wang and Jia Gou., Accumulation of the attenuation operator t* along wavefront evolution in heterogeneous absorbing media using a modified fast marching method, Computers & Geoscience,******************************.

Requirements

| Requirement  | Version   |
| ------------ | --------- |
| Python       | ≥3        |
| SCons        | latest    |
| g++          | ≥5        |
| OpenMP       | supported |
| C++ standard | C++11     |


The main contents

Within the folder 'Src_Read', there are several subroutines, we give a detail descriptions in follows.

'vel_model_load.cpp' - load the 3D velocity model in specific format.
'att_model_load.cpp' - load the 3D attenuation model (the quality factor Q) in the same format with 'vel_model_load.cpp'.
'source_load.cpp' - load the source information (For example, the location of the earthquake, or other sources).
'receiver_load.cpp' - load the receiver information (For example, the location of the station, or other receivers).
'raytracing.cpp' - tracing the ray path along the negative gradient direction of the traveltime field.
'read_parafile.cpp' - read these related parameters during the calculation.
'eikonal.cpp' - the real fast marching algorithm within and the way to calculate t* along t. 
'main.cpp' - combine the above mentioned subroutines.

Preparation

'Par_file' - the parameter file
'velocity3d' - the 3D velocity model file
'attenuation3d' - the 3D attenuation (Qp) model file
'sources' - the source location file
'receivers' - the receivers location file

Further details:

while using this package, three parameters in main.cpp need to be changed based on your object:
'usesecond = false' - determine if use the second-order derivatives during the calculation, could be true or false.
'usecross = false' - determine if use the cross derivatives, keep it set as false.
output_ts = true' - determine if output the final traveltime field and the attenuation operator t* field.







