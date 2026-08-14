# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.
#
# Regression test: runs a small, fast, deterministic simulation (a scaled-down
# version of python/examples/collapse.py - a few particles, 3 frames, single
# thread for determinism) entirely through the Python bindings and compares the
# resulting particle positions/velocities against a frozen "golden" snapshot
# checked into python/tests/golden/collapse_mini_{2,3}d.npz - one file per
# dimensionality, since THREEDIM (src/tools.hpp) changes what this scenario
# even means (different Lz/gravity/box shape) and switching it shouldn't
# invalidate the golden data for whichever dimension you're not currently in.
#
# Why frozen values instead of building examples/collapse.cpp and diff-ing live
# output each run: the golden values were originally verified byte-for-byte
# identical against a standalone C++ build of the same scenario (same compiled
# simulation/*.cpp code, just linked into a different target) - see the
# conversation this was introduced in for that one-off cross-language check.
# Once the *bindings* are known to be a straightforward, order-correct
# passthrough (which test_object_construction_argument_order.py checks
# directly), there is no ongoing risk of the Python and C++ paths silently
# diverging from each other - only of a future code change (in the bindings OR
# in simulation/*.cpp, which both targets share) silently changing behavior.
# A frozen golden snapshot catches exactly that, without needing to compile and
# run a second C++ binary on every test invocation.
#
# If you make an intentional change to the simulation algorithm (not just the
# bindings), this test will correctly start failing - regenerate the golden
# file by re-running this scenario and re-saving. The bottom of this file is
# runnable directly for that:
#   PYTHONPATH=build/python python3 python/tests/test_regression.py

from pathlib import Path

import numpy as np
import pytest

import matter

DIM = matter.Simulation().dim
GOLDEN_PATH = Path(__file__).resolve().parent / "golden" / f"collapse_mini_{DIM}d.npz"


def vec(x=0.0, y=0.0, z=0.0):
    return [x, y, z][:DIM]


def build_mini_collapse_sim():
    """A scaled-down version of python/examples/collapse.py: same setup shape
    (Hencky elasticity, DPVisc plasticity, a bottom+left ObjectPlate boundary,
    one ObjectBox obstacle) but a coarser sampling radius, fewer frames, and a
    single thread, so it runs in milliseconds and produces exactly the same
    particle state every time."""
    sim = matter.Simulation()
    sim.initialize(False)  # no disk I/O
    sim.n_threads = 1      # deterministic
    sim.end_frame = 3
    sim.fps = 10
    sim.cfl = 0.5
    sim.flip_ratio = -0.95
    sim.calculate_energy = True
    sim.use_sparse = False

    sim.elastic_model = matter.ElasticModel.Hencky
    sim.E = 1e6
    sim.nu = 0.3
    sim.rho = 1000

    sim.gravity = vec(0.0, -9.81)

    sim.Lx = 0.3
    sim.Ly = 0.3
    if DIM == 3:
        sim.Lz = 0.1
    matter.sample_particles(sim, 0.02, seed=42)

    positions = sim.particles.get_positions()
    positions[:, 0] -= 0.5 * sim.Lx
    positions[:, 1] += 0.5 * sim.dx
    sim.particles.set_positions(positions)

    sim.add_plate(matter.ObjectPlate(0, matter.PlateType.bottom, matter.BC.NoSlip))
    sim.add_plate(matter.ObjectPlate(0, matter.PlateType.left, matter.BC.SlipFree))

    box_L = vec(0.06, 0.06, 0.06)
    box_c = vec(0.3, 0.03, 0.0)
    box_v = vec(0.0, 0.0, 0.0)
    sim.add_object(matter.ObjectBox(matter.BC.SlipFree, 0.2, "object_box", True, box_L, box_c, box_v))

    sim.plastic_model = matter.PlasticModel.DPVisc
    sim.use_pradhana = True
    sim.q_prefac = 1.0 / np.sqrt(2.0)
    sim.M = np.tan(30 * np.pi / 180.0)
    sim.q_cohesion = 0
    sim.visc_exponent = 1
    sim.visc_time = 0

    return sim


@pytest.mark.skipif(not GOLDEN_PATH.exists(),
                     reason=f"no {DIM}D golden file - see this file's docstring to generate one")
def test_collapse_mini_matches_golden():
    sim = build_mini_collapse_sim()
    sim.simulate()
    assert sim.exit == 0

    golden = np.load(GOLDEN_PATH)
    assert int(golden["dim"][0]) == DIM  # sanity check against a mislabeled/miscopied golden file
    assert sim.Np == int(golden["Np"][0])
    assert sim.dx == pytest.approx(float(golden["dx"][0]), rel=1e-8)

    positions = sim.particles.get_positions()
    velocities = sim.particles.get_velocities()
    assert positions == pytest.approx(golden["positions"], rel=1e-8, abs=1e-12)
    assert velocities == pytest.approx(golden["velocities"], rel=1e-8, abs=1e-12)


if __name__ == "__main__":
    sim = build_mini_collapse_sim()
    sim.simulate()
    assert sim.exit == 0
    GOLDEN_PATH.parent.mkdir(exist_ok=True)
    np.savez(
        GOLDEN_PATH,
        dim=np.array([DIM]),
        Np=np.array([sim.Np]),
        dx=np.array([sim.dx]),
        positions=sim.particles.get_positions(),
        velocities=sim.particles.get_velocities(),
    )
    print(f"Wrote {GOLDEN_PATH} (Np={sim.Np}, dim={DIM})")
