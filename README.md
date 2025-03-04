# Matter

[![Build and Test](https://github.com/larsblatny/matter/actions/workflows/cmake-single-platform.yml/badge.svg)](https://github.com/larsblatny/matter/actions/workflows/cmake-single-platform.yml)

_Matter_ is an open-source C++ implementation of the Material Point Method (MPM) with elasto-viscoplastic rheologies, specifically designed to model the mechanics and flow of granular matter. However, its usage extends to simulate a variety of different _matter_ undergoing small and large deformations.    

Developed across different sides of the Swiss Alps, this software is designed to be lightweight, easy to install and use, with few dependencies, without significantly compromising speed and performance. Parallelized on shared memory with OpenMP, millions of material points/particles can be simulated on desktops and laptops.

![logo](https://larsblatny.github.io/images/matter.png)

## Examples
With _Matter_, you can simulate granular flow on various simple and complex terrains, as demonstrated below.

<p align="middle">
  <img src="https://larsblatny.github.io/images/double_bump_cropped_medium.gif" height="250" /> 
  <img src="https://larsblatny.github.io/images/mountain_cropped_medium.gif" height="250" />
</p>

<p align="middle">
  <img src="https://larsblatny.github.io/images/collapse_fps50_cropped_medium.gif" height="262" />
  <img src="https://larsblatny.github.io/images/silo_fps30_10X_cropped_medium.gif" height="262" /> 
</p>

## License and attribution
Matter is an open-source software licensed under _GNU General Public License v3.0_ (see LICENSE file).
If you are interested in using Matter in commercial products or services, please do not hesitate to contact [Lars Blatny](https://larsblatny.github.io/) (lars.blatny [at] slf.ch).

If you use Matter in your research, please cite the scientific works where this code has been used:   
* Blatny, L., Gray, J.M.N.T. and Gaume, J. (2024) _A critical state μ(I)-rheology model for cohesive granular flows_, Journal of Fluid Mechanics, 997, p. A67. [doi:10.1017/jfm.2024.643](https://doi.org/10.1017/jfm.2024.643)


## Features

* Both **2D** (plane strain) and **3D**


* Supports Update Stress Last **(USL)** and Modified Update Stress Last **(MUSL)** updated Lagrangian **symplectic** (explicit) MPM


* Particle <-> grid interplation with **quadratic** and **cubic B-splines** with various transfer schemes:
    *  Linearly weighted combination of **PIC** and **FLIP** (Stomakhin et al., 2013)
    *  **APIC** (Jiang et al., 2015)
    *  **AFLIP** (Fei et al., 2021)


* Adaptive background **grid extension/reduction** to particle domain


* **Finite strain** elastoplastic framework for constitutive models
  

* The following **(visco)plastic models** are currently available and compatible with Hencky's elasticity model:
    * Von Mises
    * Drucker-Prager
    * Modified Cam-Clay
    * Perzyna-Von Mises
    * Perzyna-Drucker-Prager
    * Perzyna-Modified Cam-Clay
    * mu(I)-Drucker-Prager
    * mu(I)-Modified Cam-Clay
    * Implementing other (visco)plastic models is easy due to the general framework of the code


* Supports **analytical objects** and **complex terrains**


* Supports no-slip as well as frictional **boundary conditions** which can be supplied a **Coulomb friction parameter**
    * No-slip (`NOSLIP`)
    * Frictional slip free to separate (`SLIPFREE`)
    * Frictional slip constrained to boundary (`SLIPSTICK`)   


* **Material-induced boundary friction (MIBF)**: A potentially variable internal friction parameter of a particular plastic model is used as the Coulomb friction for terrain-material interaction

* **Parallelized** on shared memory with [OpenMP](https://www.openmp.org/)

* Initial particle positions can be sampled using the **Poisson disk sampling** scheme by R. Bridson, ACM SIGGRAPH 2007, here based on the implementation by [Tommy Hinks](https://github.com/thinks/poisson-disk-sampling)


## Get started

### Installing dependencies

The required non-header-only dependencies are **[CMake](https://cmake.org/)**, **[OpenMP](https://www.openmp.org/)**. The standard C++ linear algebra template library **[Eigen](https://eigen.tuxfamily.org/)** is required. The option `-DUSE_VDB=ON` also requires **[OpenVDB](https://www.openvdb.org/)**, however, this can be turned off if only analytic objects are used. 

##### Linux
Install all dependencies with this simple command:        
`sudo apt-get install -y cmake libeigen3-dev libopenvdb-dev libtbb-dev libboost-all-dev libilmbase-dev libopenexr-dev`    
where only `cmake` and `libeigen3-dev` are required. 

##### MacOS
You can download/install the required dependencies with [Homebrew](https://brew.sh/) through   
`brew install cmake eigen libomp`    
and OpenVDB with (although this may require other dependecies, see https://formulae.brew.sh/formula/openvdb)    
`brew install openvdb`    

However, we recommend instead using [VSCode Dev Containers](https://code.visualstudio.com/docs/devcontainers/tutorial). 
Download [VSCode](https://code.visualstudio.com/) and [Docker](https://www.docker.com/products/docker-desktop/). In VSCode, install the "Dev Containers" extension. When you open the Matter repository in VSCode, you will be asked if you want to run it in a container. When you accept, a container will be set up and you are ready to simulate!


### Build and run the code

1. Set up your simulation parameters and initial state in `mpm.cpp`.  In `tools.hpp`, the user can specify the dimension of the simulation (default: 2D) and order of interpolation (default: quadratic). 

2. Create a build directory:     
   `mkdir build`

3. Enter the build directory:    
   `cd build`
   
4. Specify CMake options:     
   `cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VDB=ON ..`     
   or    
   `cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VDB=OFF ..`       

5. Compile (NB: the number of threads for the _simulation_ is specified in `mpm.cpp`)      
   `make -j <number of cores for compilation>` 

6. Run the executable:      
   `./src/mpm`

### Example of setup file

Here is a minimal setup file `mpm.cpp` for a simple granular collapse.    

```cpp
#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"
#include "objects/object_plate.hpp"

int main(){

    Simulation sim;

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "collapse");

    sim.end_frame = 20;     // last frame to simulate
    sim.fps = 10;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel

    sim.Lx = 1;
    sim.Ly = 1;
    sim.Lz = 1; // ONLY IF 3D, OTHERWISE REMOVE LINE
    SampleParticles(sim, /*sampling radius*/ 0.01);

    ObjectPlate ground = ObjectPlate(/*position*/ 0, /*type*/ bottom, /*boundary cond.*/ SLIPFREE, /*friction*/ 0.5);  
    sim.plates.push_back(ground);

    sim.rho = 1000;         // Density (kg/m3)
    sim.gravity[1] = -9.81; // Gravity

    sim.E = 1e6;    // Young's modulus (Pa)
    sim.nu = 0.3;   // Poisson's ratio (-)

    sim.plastic_model = DP; // Drucker-Prager yield
    sim.dp_slope = 0.5;     // Internal friction
    sim.dp_cohesion = 0;    // Cohesion

    sim.simulate();
    return 0;
}
```

In the `examples` folder, other examples will be archived (to use one of these examples, simply copy it into the `src` folder and rename it `mpm.cpp`).

### Objects and terrains

Rigid objects and terrains (boundaries) are either   
* formulated analytically as level sets (signed distance functions)   
* or imported as `.vdb` level sets files using [OpenVDB](https://www.openvdb.org/)   

Analytical objects can be specified as a derived class from the general `ObjectGeneral` class. An example of this is `ObjectBump` which provides the terrain of a smooth bump used in the flow experiments in [Viroulet et al. (2017)](https://doi.org/10.1017/jfm.2017.41). For the very common case of an axis-aligned plate, an `ObjectPlate` class has been made separate from `ObjectGeneral` class for efficiency and convenience. In `ObjectPlate`, you can also assign a speed to the plate, as well as controls on the time-evolution of the speed. Any plate must either a `top`, `bottom`, `front`, `back`, `left` or `right`. 

A terrain/object from a `.vdb` is stored in an instance of the `ObjectVdb` class which is derived from `ObjectGeneral`.
Examples of `.vdb`-files are found in the folder `levelsets`.

Note that all `ObjectGeneral` instances must be added to the vector `objects` and `ObjectPlate` instances are added to the vector `plates`.


### Saving simulation data

The directory to save the output data is specified by the user in the setup file `mpm.cpp`.
Particle data is saved as **binary PLY-files** (using [tinyply](https://github.com/ddiakopoulos/tinyply)) with the format (`particles_fX.ply`) where X represents the frame number (from 0 to `end_frame` as specified by the user).

We recommend [SideFX's Houdini](https://www.sidefx.com) for visualization the particle data. In the file `visualize.hipnc` in the `postprocess` directory, we show how to make a simple visualization of the data. Some simple post-processing can also be done directly in Houdini. PLY files can also be easily read by Python. This is shown in the file `load_ply.py` which can be found in the `postprocess` directory.

Optionally, the grid data can be saved if `save_grid = true`. Then, the grid data is saved as `grid_fX.ply`. By default, the grid data is not saved.


### Key parameters and options
This is a non-exhaustive list of parameters and options (of the `Simulation` class) to be specified in the input file `mpm.cpp`. See `simulation.hpp` for the complete list, and take advantage of the current `mpm.cpp` example file. Other example files are found in the `examples` directory.

| Parameter  | Default value  | Description  |
| ----       |    ----        |          ---    |
| `end_frame`  | 1   | Last frame to simulate   
| `fps`        | 1.0   |  Frames per second
| `n_threads`  | 1    |  Number of threads in parallel
| `save_grid`  | false | Save grid data to file
| `cfl`        | 0.5     | CFL constant, typically around 0.5
| `cfl_elastic`| 0.5     | CFL-like constant for elastic wave speed, typically around 0.5
| `flip_ratio` | -0.95  | (A)PIC-(A)FLIP ratio in the range [-1,1]. Positive numbers [0,1]: PIC/FLIP where 1 equals pure FLIP. Negative numbers [-1,0): APIC/AFLIP where -1 equals pure AFLIP.
| `reduce_verbose` | false | Reduce writing to screen
| `gravity`    | (0,0,(0)) | Gravitational acceleration vector. If `gravity_special = true` one may code the gravity evolution in `update_dt.cpp`.
| `rho`                   |  1000           | Density (kg/m3)
| `musl`                 | false          | Use MUSL instead of USL
| `use_mibf` | false           | Use Material-Induced Boundary Friction (MIBF), only relevant for certain plasticity models
| `use_von_mises_q`      | false          | If `true` define the "equivalent shear stress" q as the von Mises equivalent stress q = sqrt(3/2 s:s), otherwise q = sqrt(1/2 s:s). 
| `pbc`        | false | Use periodic boundary conditions, see `pbc.cpp`. If `pbc_special = true` one may code the periodicity in `position_update.cpp`
| `Lx`, `Ly` and `Lz`    | 1.0            | The material sample space used in `SampleParticles(...)`
| `grid_reference_point` | - | Optionally provide a point to be considered in the initial adaptive grid creation, otherwise it only considers the particle domain
| `elastic_model`        | Hencky | Elastic model. Note that Hencky's model must be used when combined with a plastic model. 
| `plastic_model`        | NoPlasticity   | Plastic model. Parameters are set according to the model used, see below.
| `E`                    | 1e6   | The 3D Young's modulus (Pa)
| `nu`                   | 0.3  | The 3D Poisson's ratio (-)


Here is a list of the various plastic models and their parameters:

| Model                                | Name                  | Parameters                | Default value   |
|  ----                                | ----     |    ----        |          ---    |
| Von Mises                            | `VM`   | `q_max`        | 100.0           |
| Drucker-Prager                       | `DP`   | `dp_slope`     | 1.0             |
|                                      |          | `dp_cohesion`  | 0.0             |
| Drucker-Prager with strain-softening | `DPSoft` | `dp_slope`   | 1.0             |
|                                      |            | `dp_cohesion`  | 0.0           |
|                                      |            | `xi`           | 0.0           |  
|                                      |            | `use_pradhana` | true          |  
| Modified Cam-Clay | `MCC`  | `beta`         | 0.0             |
|                   |          | `p0`           | 1000.0          |
|                   |          | `xi`           | 0.0             |
|                   |          | `M`            | 1.0             |      
| Perzyna-Von Mises | `VMVisc` | `q_max` | 100.0         |
| |                     | `q_min`        | 100.0           |
| |                     | `p_min`        | -1.0e20         |
| |                     | `xi`           | 0.0             |  
| |                     | `perzyna_exp`  | 1.0             |      
| |                     | `perzyna_visc` | 0.0             |  
| Perzyna-Drucker-Prager | `DPVisc` | `dp_slope`  | 1.0  |
| |                     | `dp_cohesion`  | 0.0             |
| |                     | `use_pradhana` | true            |
| |                     | `perzyna_exp`  | 1.0             |      
| |                     | `perzyna_visc` | 0.0             |
| Perzyna-Modified Cam-Clay | `MCCVisc`  | `beta`  | 0.0 |
| |                     | `p0`           | 1000.0          |
| |                     | `xi`           | 0.0             |
| |                     | `M`            | 1.0             |  
| |                     | `perzyna_exp`  | 1.0             |      
| |                     | `perzyna_visc` | 0.0             |
| $\mu(I)$-rheology     | `DPMui`  | `dp_cohesion` | 0.0 |
| |                     | `use_pradhana` | true            |
| |                     | `rho_s`        | 1.0             |      
| |                     | `grain_diameter`| 0.001          |
| |                     | `I_ref`        | 0.279           |      
| |                     | `mu_1`         | 0.382           |      
| |                     | `mu_2`         | 0.644           |      
| Critical state $\mu(I)$-rheology | `MCCMui`  | `beta` | 0.0  |
| |                     | `p0`           | 1000.0          |
| |                     | `xi`           | 0.0             |
| |                     | `rho_s`        | 1000.0          |      
| |                     | `grain_diameter`| 0.001          |
| |                     | `I_ref`        | 0.279           |      
| |                     | `mu_1`         | 0.382           |      
| |                     | `mu_2`         | 0.644           |

In the MCC-based models, one must also choose a corresponding hardening law: 
* exponential explicit law `ExpoExpl`    
* hyperbolic sine explicit law `SinhExpl`     
* exponential implicit law `ExpoImpl`     
* hyperbolic sine implicit law `SinhImpl`     


## Limitations

* Optimized for dense (not sparse) topologies. On sparse topologies, the simulations can be rather slow. This follows from the chosen grid data format, which may be replaced in the future. However, the current choice keeps the code simple and organized.

* Supports only single-materials, however, one can easily extend this to multi-material problems. E.g., one can create particle quantities for the relevant material parameters (see `data_structures.hpp`) which can then be used in the material models (see, e.g., `plasticity.cpp`)

## Troubleshooting

* On Mac, especially if OpenMP was installed through Homebrew, you might need to specify certiain paths for CMake, e.g.,       
`cmake -DCMAKE_BUILD_TYPE=Release -DOpenMP_CXX_FLAG="-Xclang -fopenmp" -DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_CXX_LIB_NAMES=libomp -DOpenMP_C_FLAG="-Xclang -fopenmp" -DOpenMP_C_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_C_LIB_NAMES=libomp -DOpenMP_libomp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib ..`

* Eigen error? Remember to specify vectors with 3 elements for 3D problems and 2 elements for 2D problems. The dimension of the problem is chosen as a global variable in `tools.hpp`.

## Help? Want to contribute?

Please contact [Lars Blatny](https://larsblatny.github.io/) (lars.blatny [at] slf.ch)
