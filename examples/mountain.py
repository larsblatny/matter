# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

import numpy as np
import matter

sim = matter.Simulation()

sim.initialize(True, "output/", "mountain")

sim.save_grid = False
sim.end_frame = 100      # last frame to simulate
sim.fps = 5               # frames per second
sim.n_threads = 8         # number of threads in parallel
sim.cfl = 0.5              # CFL constant, typically around 0.5
sim.flip_ratio = -0.95    # (A)PIC-(A)FLIP ratio in [-1,1].
sim.reduce_verbose = True

sim.use_sparse = True

# INITIALIZE ELASTICITY
sim.elastic_model = matter.ElasticModel.Hencky
sim.E = 5e5      # Young's modulus (Pa)
sim.nu = 0.3     # Poisson's ratio (-)
sim.rho = 2000   # Density (kg/m3)

# GRAVITY ANGLE [default: gravity is 0]
theta_deg = 30.0   # angle in degrees of gravity vector
theta = theta_deg * np.pi / 180
sim.gravity = [9.81 * np.sin(theta), -9.81 * np.cos(theta), 0]

# INITIAL PARTICLE POSITIONS
release_path = "levelsets/mountain_release.vdb"
release = matter.ObjectVdb(release_path)
matter.sample_particles_from_vdb(sim, release, 0.18)             # h = 0.431492478463 m
# matter.sample_particles_from_vdb(sim, release, 0.1668627)      # h = 0.400046583858 m
# matter.sample_particles_from_vdb(sim, release, 0.06663595279)  # h = 0.159996934438 m

terrain_path = "levelsets/mountain_terrain_20mb.vdb"
sim.add_object(matter.ObjectVdb(terrain_path, matter.BC.SlipFree, 0.35, "Terrain"))

sim.plastic_model = matter.PlasticModel.DP   # Perzyna model with Drucker-Prager yield surface

sim.q_prefac = 1.0 / np.sqrt(2.0)   # [default: sqrt(1/2)] Prefactor in def. of q

sim.M = np.tan(30 * np.pi / 180.0)   # Internal friction
sim.q_cohesion = 0   # Yield surface's intersection of q-axis (Pa), 0 is cohesionless

sim.simulate()
