# MFMM_tstar
MFMM_tstar

What is MFMM_tstar ?

tstar_FMM is a software for calculating attenuation operator t* directly by the fast marching method without ray tracing.

*** Please cite the software as:*** Dongdong wang and Jia Gou., Accumulation of the attenuation operator t* along wavefront evolution in heterogeneous absorbing media using a modified fast marching method, Computers & Geoscience,******************************.

How to get set up?

what the package includes ?

Within the folder 'Src_Read', there are several subroutines, we give a detail descriptions in follows.

'vel_model_load.cpp' - load the 3D velocity model in specific format.
'att_model_load.cpp' - load the 3D attenuation model (the quality factor Q) in the same format with 'vel_model_load.cpp'.
