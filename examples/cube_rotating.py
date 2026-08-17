# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

import numpy as np
import matter

DIM = matter.Simulation().dim
def vec(x=0.0, y=0.0, z=0.0):
    return [x, y, z][:DIM]

sim = matter.Simulation()

sim.initialize(True, "output/", "cube_rotating")
sim.save_grid = True
sim.reduce_verbose = True

sim.end_frame = 20
sim.fps = 1

sim.gravity = vec()

sim.cfl = 0.5
sim.flip_ratio = -1
sim.n_threads = 8

sim.Lx = 1
sim.Ly = 1
if DIM == 3:
    sim.Lz = 0.05
matter.sample_particles(sim, 0.01)

positions = sim.particles.get_positions()
positions[:, 0] -= 0.5 * sim.Lx
positions[:, 1] -= 0.5 * sim.Ly

vx = -1.0 * positions[:, 1] + 0.5
vy = 1.0 * positions[:, 0] + 0.5
velocity_cols = [vx, vy] + ([np.zeros_like(vx)] if DIM == 3 else [])
velocities = np.column_stack(velocity_cols)

total_energy_init = 0.5 * np.sum(vx * vx + vy * vy)  # per unit mass

sim.particles.set_positions(positions)
sim.particles.set_velocities(velocities)

# Elasticity
sim.elastic_model = matter.ElasticModel.Hencky
sim.plastic_model = matter.PlasticModel.NoPlasticity
sim.E = 1e6      # Young's modulus (Pa)
sim.nu = 0.3     # Poisson's ratio (-)
sim.rho = 1550   # Density (kg/m3)

sim.simulate()

final_velocities = sim.particles.get_velocities()
total_energy_last = 0.5 * np.sum(np.sum(final_velocities ** 2, axis=1))  # per unit mass

print("E_init =", total_energy_init)
print("E_last =", total_energy_last)

rel_diff = (total_energy_init - total_energy_last) / total_energy_init
print("rel_diff =", rel_diff)
