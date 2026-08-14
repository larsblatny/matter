# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.
#
# Tests for the pybind11 bindings themselves (python/bindings.cpp) - NOT for MPM
# physics, which is already covered by the C++ gtest suite in tests/tests.cpp
# (identical compiled simulation code runs either way, so re-testing physics here
# would just duplicate that suite through an extra layer). These tests target
# things that are unique to the Python binding layer and that the C++ tests
# structurally cannot catch: enum value mappings, constructor argument order and
# defaults, object copy/dispatch semantics through add_object/add_plate, numpy
# <-> Eigen conversion behavior, and error propagation from C++ into Python.
#
# Run with (from the repo root, after building):
#   PYTHONPATH=build/python python3 -m pytest python/tests -v

from pathlib import Path

import numpy as np
import pytest

import matter

DIM = matter.Simulation().dim
LEVELSETS_DIR = Path(__file__).resolve().parents[2] / "levelsets"


def vec(x=0.0, y=0.0, z=0.0):
    return [x, y, z][:DIM]


# --------------------------------------------------------------------------- #
# Module / construction
# --------------------------------------------------------------------------- #

def test_module_imports_and_constructs():
    sim = matter.Simulation()
    assert sim.dim in (2, 3)
    assert sim.exit == 0


def test_dim_is_readonly():
    sim = matter.Simulation()
    with pytest.raises(AttributeError):
        sim.dim = 3


# --------------------------------------------------------------------------- #
# Enums: distinct values, and round-trip through bound fields
# --------------------------------------------------------------------------- #

def test_elastic_model_enum_roundtrips():
    sim = matter.Simulation()
    for value in (matter.ElasticModel.Hencky, matter.ElasticModel.NeoHookean):
        sim.elastic_model = value
        assert sim.elastic_model == value


def test_plastic_model_enum_values_are_distinct_and_roundtrip():
    sim = matter.Simulation()
    values = [
        matter.PlasticModel.NoPlasticity, matter.PlasticModel.VM, matter.PlasticModel.DP,
        matter.PlasticModel.DPSoft, matter.PlasticModel.MCC, matter.PlasticModel.VMVisc,
        matter.PlasticModel.DPVisc, matter.PlasticModel.MCCVisc, matter.PlasticModel.DPMui,
        matter.PlasticModel.MCCMui,
    ]
    assert len(set(values)) == len(values)  # no accidental aliasing between enumerators
    for value in values:
        sim.plastic_model = value
        assert sim.plastic_model == value


def test_hardening_law_enum_values_are_distinct():
    values = [matter.HardeningLaw.NoHard, matter.HardeningLaw.ExpoExpl,
              matter.HardeningLaw.ExpoImpl, matter.HardeningLaw.SinhExpl,
              matter.HardeningLaw.SinhImpl]
    assert len(set(values)) == len(values)


def test_bc_enum_roundtrips_through_object_construction():
    for value in (matter.BC.NoSlip, matter.BC.SlipStick, matter.BC.SlipFree):
        obj = matter.ObjectGround(value, 0.1)
        assert obj.bc == value


def test_plate_type_enum_roundtrips_through_object_plate():
    for value in (matter.PlateType.top, matter.PlateType.bottom, matter.PlateType.left,
                  matter.PlateType.right, matter.PlateType.front, matter.PlateType.back):
        plate = matter.ObjectPlate(0.0, value)
        assert plate.plate_type == value


# --------------------------------------------------------------------------- #
# Simulation scalar field readwrite
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("field,value", [
    ("n_threads", 4), ("end_frame", 7), ("fps", 12.5), ("cfl", 0.3),
    ("flip_ratio", -0.5), ("rho", 2100.0), ("E", 2e6), ("nu", 0.25),
    ("M", 0.7), ("q_cohesion", 10.0), ("use_sparse", True), ("save_grid", True),
    ("calculate_energy", True), ("use_pradhana", False),
])
def test_scalar_fields_roundtrip(field, value):
    sim = matter.Simulation()
    setattr(sim, field, value)
    assert getattr(sim, field) == value


# --------------------------------------------------------------------------- #
# Eigen vector field semantics - a real, documented gotcha
# --------------------------------------------------------------------------- #

def test_gravity_whole_vector_assignment_persists():
    sim = matter.Simulation()
    sim.gravity = vec(1.0, -2.0)
    assert list(sim.gravity) == pytest.approx(vec(1.0, -2.0))


def test_gravity_item_assignment_raises_instead_of_silently_failing():
    """Reading an Eigen-typed field returns a *copy* (pybind11/eigen.h default
    behavior for a plain-value member), and pybind11 marks that copy read-only,
    so `sim.gravity[0] = x` raises immediately rather than silently discarding
    the write. Good: it fails loudly instead of losing data. Always assign the
    whole vector instead (`sim.gravity = [gx, gy]`) to actually mutate the field.
    If this test starts failing, the getter's semantics changed (e.g. now
    returns a writable view) and the example scripts' comments need revisiting."""
    sim = matter.Simulation()
    sim.gravity = vec(0.0, 0.0)
    with pytest.raises(ValueError):
        sim.gravity[0] = 99.0


# --------------------------------------------------------------------------- #
# Particles
# --------------------------------------------------------------------------- #

def test_particles_default_construction():
    p = matter.Particles(5)
    assert len(p.x) == 5
    assert len(p.F) == 5
    assert list(p.x[0]) == pytest.approx(vec())
    assert np.array(p.F[0]) == pytest.approx(np.eye(DIM))
    assert p.eps_pl_dev[0] == 0.0


def test_particles_list_of_array_field_roundtrip():
    p = matter.Particles(3)
    new_x = [vec(1, 2), vec(3, 4), vec(5, 6)]
    p.x = new_x
    for got, expected in zip(p.x, new_x):
        assert list(got) == pytest.approx(expected)


def test_particles_get_set_positions_bulk_roundtrip():
    p = matter.Particles(4)
    arr = np.arange(4 * DIM, dtype=float).reshape(4, DIM)
    p.set_positions(arr)
    got = p.get_positions()
    assert got.shape == (4, DIM)
    assert got == pytest.approx(arr)


def test_particles_set_positions_wrong_row_count_raises():
    p = matter.Particles(4)
    with pytest.raises(RuntimeError):
        p.set_positions(np.zeros((3, DIM)))


# --------------------------------------------------------------------------- #
# sample_particles / sample_particles_multi
# --------------------------------------------------------------------------- #

def _set_domain(sim, Lx=1.0, Ly=1.0, Lz=1.0):
    sim.Lx = Lx
    sim.Ly = Ly
    if DIM == 3:
        sim.Lz = Lz


def test_sample_particles_populates_fields_within_domain():
    sim = matter.Simulation()
    _set_domain(sim)
    matter.sample_particles(sim, 0.05, seed=1)
    assert sim.Np > 0
    assert len(sim.particles.x) == sim.Np
    assert sim.dx > 0
    assert sim.particle_mass > 0
    upper = vec(sim.Lx, sim.Ly, getattr(sim, "Lz", 0.0))
    for pos in sim.particles.x:
        for d in range(DIM):
            assert -1e-9 <= pos[d] <= upper[d] + 1e-9


def test_sample_particles_is_deterministic_given_a_seed():
    sim1, sim2 = matter.Simulation(), matter.Simulation()
    _set_domain(sim1)
    _set_domain(sim2)
    matter.sample_particles(sim1, 0.05, seed=7)
    matter.sample_particles(sim2, 0.05, seed=7)
    assert sim1.Np == sim2.Np
    for a, b in zip(sim1.particles.x, sim2.particles.x):
        assert list(a) == pytest.approx(list(b))


def test_sample_particles_multi_populates_start_indices():
    sim = matter.Simulation()
    origins = [vec(0, 0), vec(2, 0)]
    sizes = [vec(1, 1), vec(1, 1)]
    matter.sample_particles_multi(sim, origins, sizes, 0.05, seed=3)
    assert sim.Np > 0
    assert len(sim.sampling_start_idx) == len(origins) + 1


# --------------------------------------------------------------------------- #
# ObjectPlate
# --------------------------------------------------------------------------- #

def test_object_plate_defaults():
    plate = matter.ObjectPlate(0.5, matter.PlateType.bottom)
    assert plate.pos_object == pytest.approx(0.5)
    assert plate.bc == matter.BC.NoSlip
    assert plate.vmin_factor == pytest.approx(1.0)
    assert plate.force_calc is False


def test_add_plate_stores_a_copy_not_an_alias():
    sim = matter.Simulation()
    plate = matter.ObjectPlate(0.1, matter.PlateType.left, matter.BC.SlipFree, 0.2)
    sim.add_plate(plate)
    assert len(sim.plates) == 1
    assert sim.plates[0].pos_object == pytest.approx(0.1)
    assert sim.plates[0].bc == matter.BC.SlipFree
    sim.plates[0].pos_object = 999.0
    assert plate.pos_object == pytest.approx(0.1)  # original untouched


# --------------------------------------------------------------------------- #
# ObjectGeneral hierarchy / add_object polymorphic dispatch
# --------------------------------------------------------------------------- #

def _make_one_of_each_object():
    objs = [
        matter.ObjectBox(matter.BC.SlipFree, 0.2, "box", False, vec(1, 1), vec(0, 0), vec(0, 0)),
        matter.ObjectBoxRotated(matter.BC.SlipFree, 0.2, "boxrot", False,
                                 vec(1, 1), vec(0, 0), vec(0, 0), 0.1),
        matter.ObjectBump(matter.BC.NoSlip, 0.1),
        matter.ObjectCurve(matter.BC.NoSlip, 0.1),
        matter.ObjectGate(matter.BC.NoSlip, 0.1),
        matter.ObjectGround(matter.BC.NoSlip, 0.1),
        matter.ObjectGroundRotated(matter.BC.NoSlip, 0.1),
        matter.ObjectRamp(matter.BC.NoSlip, 0.1),
    ]
    if DIM == 3:
        objs.append(matter.ObjectSilo(matter.BC.NoSlip, 0.1))
    return objs


def test_add_object_dispatches_to_correct_derived_python_type():
    sim = matter.Simulation()
    objs = _make_one_of_each_object()
    for o in objs:
        sim.add_object(o)
    assert len(sim.objects) == len(objs)
    for original, stored in zip(objs, sim.objects):
        assert type(stored) is type(original)
        assert isinstance(stored, matter.ObjectGeneral)


def test_stored_object_exposes_base_and_derived_fields():
    sim = matter.Simulation()
    box = matter.ObjectBox(matter.BC.SlipFree, 0.3, "mybox", True, vec(2, 2), vec(1, 1), vec(0, 0))
    sim.add_object(box)
    stored = sim.objects[0]
    assert stored.name == "mybox"
    assert stored.friction == pytest.approx(0.3)
    assert stored.force_calc is True
    # L is a derived-only field; only reachable if pybind11 returned the real subtype
    assert list(stored.L) == pytest.approx(vec(2, 2))


def test_add_object_rejects_unregistered_type():
    sim = matter.Simulation()
    with pytest.raises(TypeError):
        sim.add_object("not an object")


def test_add_object_stores_a_copy_not_an_alias():
    sim = matter.Simulation()
    box = matter.ObjectBox(matter.BC.SlipFree, 0.2, "box")
    sim.add_object(box)
    sim.objects[0].friction = 0.9
    assert box.friction == pytest.approx(0.2)  # original untouched


@pytest.mark.skipif(not hasattr(matter, "ObjectVdb"), reason="built without USE_VDB")
def test_object_vdb_construction_and_add_object():
    vdb_path = str(LEVELSETS_DIR / "box.vdb")
    obj = matter.ObjectVdb(vdb_path, matter.BC.NoSlip, 0.1, "box_vdb")
    sim = matter.Simulation()
    sim.add_object(obj)
    stored = sim.objects[-1]
    assert type(stored) is matter.ObjectVdb
    assert stored.name == "box_vdb"
    min_bbox, max_bbox = stored.bounds()
    assert len(min_bbox) == DIM
    assert len(max_bbox) == DIM


# --------------------------------------------------------------------------- #
# simulate() error propagation
# --------------------------------------------------------------------------- #

def test_simulate_raises_when_exit_flag_is_set():
    sim = matter.Simulation()
    # simulate() never resets `exit` to 0, so setting it beforehand reliably
    # exercises the wrapper's post-call check without needing a real solver failure.
    sim.exit = 1
    with pytest.raises(RuntimeError):
        sim.simulate()


def test_simulate_without_initialize_is_a_noop_not_an_error():
    """Documents current behavior: forgetting to call initialize() does not raise
    from Python - the C++ simulate() just prints a message and returns early with
    exit still 0. Pinning this down means a future change to either side is a
    deliberate decision, not a silent surprise."""
    sim = matter.Simulation()
    sim.simulate()
    assert sim.exit == 0


# --------------------------------------------------------------------------- #
# End-to-end smoke test through the full binding surface (no disk I/O)
# --------------------------------------------------------------------------- #

def test_minimal_simulation_runs_end_to_end():
    sim = matter.Simulation()
    sim.initialize(False)  # save_sim=False: no directory/mpm.cpp-copy side effects
    sim.n_threads = 1
    sim.end_frame = 1
    sim.fps = 10
    _set_domain(sim, Lx=0.2, Ly=0.2, Lz=0.2)
    matter.sample_particles(sim, 0.03, seed=1)
    sim.gravity = vec(0.0, -9.81)
    sim.elastic_model = matter.ElasticModel.Hencky
    sim.E = 1e6
    sim.nu = 0.3
    sim.rho = 1000
    sim.add_plate(matter.ObjectPlate(0.0, matter.PlateType.bottom))

    sim.simulate()

    assert sim.exit == 0
    # objects/plates stay valid and readable after the run (reference_internal lifetime check)
    assert sim.plates[0].bc == matter.BC.NoSlip
