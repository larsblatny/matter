# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

import numpy as np
import matter

DIM = matter.Simulation().dim
def vec(x=0.0, y=0.0, z=0.0):
    return [x, y, z][:DIM]

sim = matter.Simulation()

sim.initialize(True, "output/", "block_incline")

sim.save_grid = True
sim.end_frame = 10       # last frame to simulate
sim.fps = 20              # frames per second
sim.n_threads = 8         # number of threads in parallel
sim.cfl = 0.5              # CFL constant, typically around 0.5
sim.flip_ratio = -0.95    # (A)PIC-(A)FLIP ratio in [-1,1].
sim.reduce_verbose = True

# INITIALIZE ELASTICITY 
sim.elastic_model = matter.ElasticModel.Hencky
sim.plastic_model = matter.PlasticModel.NoPlasticity
sim.E = 1e6      # Young's modulus (Pa)
sim.nu = 0.3     # Poisson's ratio (-)
sim.rho = 1000   # Density (kg/m3)

# INCLINE ANGLE AND FRICTION 
theta = np.radians(24)
sim.gravity = vec(9.81 * np.sin(theta), -9.81 * np.cos(theta))

# PARTICLES
block_size = vec(0.1, 0.05, 0.05 if DIM == 3 else 0.1)
sim.Lx = block_size[0]
sim.Ly = block_size[1]
if DIM == 3:
    sim.Lz = block_size[2]

k_rad = 0.005
matter.sample_particles(sim, k_rad, ppc=4 if DIM == 2 else 8)

positions = sim.particles.get_positions()
positions[:, 1] += 0.5 * sim.dx
sim.particles.set_positions(positions)

# FRICTION DATA
mu = np.tan(np.radians(15)) 
basal_friction = matter.BasalFrictionField()
if DIM == 3:
    friction_grid = np.full((2, 2), mu)  # of shape (nz, nx)
    basal_friction.set(friction_grid, x0=0.0, z0=0.0, d=max(sim.Lx, sim.Lz))
else:
    friction_line = np.full(2, mu)
    basal_friction.set(friction_line, x0=0, d=sim.Lx)
sim.basal_friction_field = basal_friction
sim.use_basal_friction_field = True

# sim.add_plate(matter.ObjectPlate(0, matter.PlateType.bottom, matter.BC.SlipFree, mu))
sim.add_object(matter.ObjectGround(matter.BC.SlipFree, mu, "ground", False, 0))

sim.simulate()

# VERIFICATION:
# under Coulomb friction: a = g*(sin(theta) - mu*cos(theta))
final_positions = sim.particles.get_positions()
final_velocities = sim.particles.get_velocities()

x_sim = final_positions[:, 0].mean()
vx_sim = final_velocities[:, 0].mean()

T_total = sim.end_frame / sim.fps
a_analytical = 9.81 * (np.sin(theta) - mu * np.cos(theta))
vx_analytical = a_analytical * T_total
dx_analytical = 0.5 * a_analytical * T_total**2

x0_mean = block_size[0] / 2   # initial mean x position (block starts at x=0)
print("Simulated slide distance:  ", x_sim - x0_mean, "m")
print("Analytical slide distance: ", dx_analytical, "m")
print("Simulated velocity:  ", vx_sim, "m/s")
print("Analytical velocity: ", vx_analytical, "m/s")
print("Relative velocity difference: ", 100 * abs(vx_sim - vx_analytical) / vx_analytical, "%")
