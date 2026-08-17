# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

import numpy as np
import matter

DIM = matter.Simulation().dim
def vec(x=0.0, y=0.0, z=0.0):
    return [x, y, z][:DIM]

sim = matter.Simulation()

sim.initialize(True, "output/", "collapse")

sim.save_grid = True
sim.end_frame = 20      # last frame to simulate
sim.fps = 10             # frames per second
sim.n_threads = 8        # number of threads in parallel
sim.cfl = 0.5             # CFL constant, typically around 0.5
sim.flip_ratio = -0.95   # (A)PIC-(A)FLIP ratio in [-1,1].

sim.calculate_energy = True   # saves energy on each particle

sim.use_sparse = False   # recommended to use True for 3D problems

# INITIALIZE ELASTICITY
sim.elastic_model = matter.ElasticModel.Hencky
sim.E = 1e6      # Young's modulus (Pa)
sim.nu = 0.3     # Poisson's ratio (-)
sim.rho = 1000   # Density (kg/m3)

# GRAVITY ANGLE [default: gravity is 0]
# Note: sim.gravity is read as a *copy*, so item assignment (sim.gravity[0] = ...) would
# raise instead of persisting - always assign the whole vector at once, which invokes the
# real setter. vec() zero-pads to match the build's dimension (2D or 3D).
theta_deg = 0   # angle in degrees of gravity vector
theta = theta_deg * np.pi / 180
sim.gravity = vec(9.81 * np.sin(theta), -9.81 * np.cos(theta))

# INITIAL PARTICLE POSITIONS
sim.Lx = 1
sim.Ly = 1
k_rad = 0.01
if DIM == 3:
    sim.Lz = 0.2
matter.sample_particles(sim, k_rad)

# OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
positions = sim.particles.get_positions()
positions[:, 0] -= 0.5 * sim.Lx
positions[:, 1] += 0.5 * sim.dx
if DIM == 3:
    positions[:, 2] -= 0.5 * sim.Lz
sim.particles.set_positions(positions)

# OPTIONAL: INITIAL PARTICLE VELOCITIES
# sim.particles.set_velocities(...)

# OBJECTS AND TERRAINS
sim.add_plate(matter.ObjectPlate(0, matter.PlateType.bottom, matter.BC.NoSlip))
sim.add_plate(matter.ObjectPlate(0, matter.PlateType.left, matter.BC.SlipFree))

box_friction = 0.2   # Friction coefficient for object
box_force = True     # Calculate force on object (can slow down simulation)
box_L = vec(0.2, 0.2, 0.2)
box_c = vec(1, 0.1, 0)
box_v = vec()        # Imposed velocity of object
sim.add_object(matter.ObjectBox(matter.BC.SlipFree, box_friction, "object_box", box_force,
                                 box_L, box_c, box_v))

# PLASTICITY
sim.plastic_model = matter.PlasticModel.DPVisc   # Viscous model with Drucker-Prager yield surface

sim.use_pradhana = True                # Suppress unwanted volume expansion in Drucker-Prager models
sim.q_prefac = 1.0 / np.sqrt(2.0)    # [default: sqrt(1/2)] Prefactor in def. of q

sim.M = np.tan(30 * np.pi / 180.0)   # Internal friction
sim.q_cohesion = 0        # Yield surface's intersection of q-axis (Pa), 0 is cohesionless
sim.visc_exponent = 1     # Exponent in viscous models
sim.visc_time = 0         # Viscous time parameter

sim.simulate()
