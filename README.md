# MFMM_tstar
MFMM_tstar

Introduction of 'MFMM_tstar' ?

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
