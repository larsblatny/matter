# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.
#
# Targets one specific, easy-to-introduce bug: the pybind11 py::init<Types...>()
# declaration in bindings.cpp listing its C++ constructor argument types in a
# different order than the real C++ constructor declares them. Since pybind11
# forwards constructor arguments purely positionally by type, a swap between
# two same-typed parameters (e.g. two TV vectors, or two doubles) still
# compiles fine and is otherwise invisible - nothing type-checks against it.
#
# Each test below constructs an object with every constructor argument set to
# a distinct, unmistakable literal value (never reusing a value - e.g. never
# leaving two TV arguments at (0,0)) and passed *positionally*, matching the
# exact order each C++ header declares. It then reads back every field and
# asserts it holds the value that was passed for that specific parameter. If
# any two same-typed parameters were ever swapped in bindings.cpp, this fails
# immediately and names exactly which field got the wrong value - no need to
# run a simulation and infer what changed from a behavioral difference.

import math
from pathlib import Path

import pytest

import matter

DIM = matter.Simulation().dim
LEVELSETS_DIR = Path(__file__).resolve().parents[2] / "levelsets"


def vec(x=0.0, y=0.0, z=0.0):
    return [x, y, z][:DIM]


def approx_vec(actual, expected):
    assert list(actual) == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# ObjectGeneral hierarchy
# --------------------------------------------------------------------------- #

def test_object_box_argument_order():
    obj = matter.ObjectBox(matter.BC.SlipFree, 0.111, "box_test", True,
                            vec(0.21, 0.22, 0.23), vec(0.31, 0.32, 0.33), vec(0.41, 0.42, 0.43))
    assert obj.bc == matter.BC.SlipFree
    assert obj.friction == pytest.approx(0.111)
    assert obj.name == "box_test"
    assert obj.force_calc is True
    approx_vec(obj.L, vec(0.21, 0.22, 0.23))
    approx_vec(obj.c, vec(0.31, 0.32, 0.33))
    approx_vec(obj.c0, vec(0.31, 0.32, 0.33))  # c0 initialized equal to c
    approx_vec(obj.v, vec(0.41, 0.42, 0.43))


def test_object_box_rotated_argument_order():
    theta = 0.81
    obj = matter.ObjectBoxRotated(matter.BC.SlipStick, 0.211, "boxrot_test", True,
                                   vec(0.51, 0.52, 0.53), vec(0.61, 0.62, 0.63),
                                   vec(0.71, 0.72, 0.73), theta)
    assert obj.bc == matter.BC.SlipStick
    assert obj.friction == pytest.approx(0.211)
    assert obj.name == "boxrot_test"
    assert obj.force_calc is True
    approx_vec(obj.L, vec(0.51, 0.52, 0.53))
    approx_vec(obj.c, vec(0.61, 0.62, 0.63))
    approx_vec(obj.c0, vec(0.61, 0.62, 0.63))
    approx_vec(obj.v, vec(0.71, 0.72, 0.73))
    assert obj.ct == pytest.approx(math.cos(theta))
    assert obj.st == pytest.approx(math.sin(theta))


def test_object_bump_argument_order():
    obj = matter.ObjectBump(matter.BC.SlipFree, 0.311, "bump_test", True)
    assert obj.bc == matter.BC.SlipFree
    assert obj.friction == pytest.approx(0.311)
    assert obj.name == "bump_test"
    assert obj.force_calc is True


def test_object_curve_argument_order():
    obj = matter.ObjectCurve(matter.BC.SlipStick, 0.411, "curve_test", True)
    assert obj.bc == matter.BC.SlipStick
    assert obj.friction == pytest.approx(0.411)
    assert obj.name == "curve_test"
    assert obj.force_calc is True


def test_object_gate_argument_order():
    obj = matter.ObjectGate(matter.BC.SlipFree, 0.511, "gate_test", True, 0.611)
    assert obj.bc == matter.BC.SlipFree
    assert obj.friction == pytest.approx(0.511)
    assert obj.name == "gate_test"
    assert obj.force_calc is True
    assert obj.height == pytest.approx(0.611)


def test_object_ground_argument_order():
    obj = matter.ObjectGround(matter.BC.SlipStick, 0.711, "ground_test", True, 0.811)
    assert obj.bc == matter.BC.SlipStick
    assert obj.friction == pytest.approx(0.711)
    assert obj.name == "ground_test"
    assert obj.force_calc is True
    assert obj.y_ground == pytest.approx(0.811)


def test_object_ground_rotated_argument_order():
    theta = 1.011
    obj = matter.ObjectGroundRotated(matter.BC.SlipFree, 0.911, "groundrot_test", True, theta)
    assert obj.bc == matter.BC.SlipFree
    assert obj.friction == pytest.approx(0.911)
    assert obj.name == "groundrot_test"
    assert obj.force_calc is True
    assert obj.ct == pytest.approx(math.cos(theta))
    assert obj.st == pytest.approx(math.sin(theta))


def test_object_ramp_argument_order():
    obj = matter.ObjectRamp(matter.BC.SlipStick, 1.111, "ramp_test", True)
    assert obj.bc == matter.BC.SlipStick
    assert obj.friction == pytest.approx(1.111)
    assert obj.name == "ramp_test"
    assert obj.force_calc is True


@pytest.mark.skipif(DIM != 3, reason="ObjectSilo is only bound for 3D builds")
def test_object_silo_argument_order():
    obj = matter.ObjectSilo(matter.BC.SlipFree, 1.211, "silo_test", True, -1.311)
    assert obj.bc == matter.BC.SlipFree
    assert obj.friction == pytest.approx(1.211)
    assert obj.name == "silo_test"
    assert obj.force_calc is True
    assert obj.cut == pytest.approx(-1.311)


@pytest.mark.skipif(not hasattr(matter, "ObjectVdb"), reason="built without USE_VDB")
def test_object_vdb_argument_order():
    vdb_path = str(LEVELSETS_DIR / "box.vdb")
    obj = matter.ObjectVdb(vdb_path, matter.BC.SlipStick, 1.411, "vdb_test", True)
    assert obj.bc == matter.BC.SlipStick
    assert obj.friction == pytest.approx(1.411)
    assert obj.name == "vdb_test"
    assert obj.force_calc is True


# --------------------------------------------------------------------------- #
# ObjectPlate - the highest-risk one: 10 (2D) / 11 (3D) positional parameters,
# several of the same type (T) right next to each other.
# --------------------------------------------------------------------------- #

def test_object_plate_argument_order():
    if DIM == 3:
        plate = matter.ObjectPlate(0.101, matter.PlateType.right, matter.BC.SlipFree, 0.202,
                                    0.303, 0.404, 0.505, 0.606, 0.657, 0.707, 0.808,
                                    "plate_test", True)
    else:
        plate = matter.ObjectPlate(0.101, matter.PlateType.right, matter.BC.SlipFree, 0.202,
                                    0.303, 0.404, 0.505, 0.606, 0.707, 0.808,
                                    "plate_test", True)

    assert plate.pos_object == pytest.approx(0.101)
    assert plate.plate_type == matter.PlateType.right
    assert plate.bc == matter.BC.SlipFree
    assert plate.friction == pytest.approx(0.202)
    assert plate.pos_lower == pytest.approx(0.303)
    assert plate.pos_upper == pytest.approx(0.404)
    assert plate.vx_object == pytest.approx(0.505)
    assert plate.vx_object_original == pytest.approx(0.505)
    assert plate.vy_object == pytest.approx(0.606)
    assert plate.vy_object_original == pytest.approx(0.606)
    if DIM == 3:
        assert plate.vz_object == pytest.approx(0.657)
        assert plate.vz_object_original == pytest.approx(0.657)
        assert plate.vmin_factor == pytest.approx(0.707)
        assert plate.load_factor == pytest.approx(0.808)
    else:
        assert plate.vmin_factor == pytest.approx(0.707)
        assert plate.load_factor == pytest.approx(0.808)
    assert plate.name == "plate_test"
    assert plate.force_calc is True
