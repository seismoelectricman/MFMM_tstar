# MFMM_tstar

## Introduction

**MFMM_tstar** is a software package for computing the seismic attenuation operator t* using a **modified fast marching method (MFMM)**.

Traditional approaches for calculating t* require explicit ray tracing followed by integration along ray paths. In contrast, **MFMM_tstar** accumulates the attenuation operator directly during the evolution of the wavefront within the fast marching framework. This approach avoids explicit ray tracing and ensures consistent accumulation of attenuation along minimum traveltime paths.

This software is designed for studies of **seismic attenuation and wave propagation in heterogeneous absorbing media**.

---

## Citation

If you use this software in your research, please cite:

Dongdong Wang and Jia Gou. *Accumulation of the attenuation operator t* along wavefront evolution in heterogeneous absorbing media using a modified fast marching method.*  Computers & Geosciences.

---

## Requirements

| Requirement | Version |
|-------------|--------|
| Python | ≥ 3 |
| SCons | Latest |
| g++ | ≥ 5 |
| OpenMP | Supported |
| C++ standard | C++11 |

---

## Compilation

Compile the code using **SCons**: 

scons (bash)

## Input Files

Before running the program, the following input files must be prepared:

| File | Description |
|-------------|--------|
| Par_file | Parameter file controlling the calculation |
| velocity3d | 3-D velocity model |
| attenuation3d | 3-D attenuation model (Qp) |
| sources | Source location file |
| receivers | Receiver location file |

## Source Code Structure

The main source files are located in the folder `Src_Read`.

- **vel_model_load.cpp**  
  Loads the 3-D velocity model.

- **att_model_load.cpp**  
  Loads the 3-D attenuation model (quality factor Q).

- **source_load.cpp**  
  Reads source information such as earthquake locations.

- **receiver_load.cpp**  
  Reads receiver information such as station locations.

- **raytracing.cpp**  
  Performs ray tracing along the negative gradient of the traveltime field.

- **read_parafile.cpp**  
  Reads parameters required for the calculation.

- **eikonal.cpp**  
  Implements the fast marching algorithm and computes the traveltime t and attenuation operator t*.

- **main.cpp**  
  Main driver program combining all modules.

## Important Parameters

When using this package, several parameters in main.cpp may need to be modified depending on your application.

- **usesecond = false**  
  Determines whether second-order derivatives are used in the fast marching update. Possible values: true or false.

- **usecross = false**  
  Determines whether cross derivatives are used. It is recommended to keep this parameter set to false.

- **output_ts = true**  
  Determines whether the traveltime field and the attenuation operator t* field are written to output files.
