// Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <memory>
#include <stdexcept>
#include <cstdint>

#include "../src/tools.hpp"
#include "../src/data_structures.hpp"
#include "../src/simulation/simulation.hpp"
#include "../src/sampling/sampling_particles.hpp"
#include "../src/sampling/regular_sampling_particles.hpp"

#include "../src/objects/object_general.hpp"
#include "../src/objects/object_plate.hpp"
#include "../src/objects/object_box.hpp"
#include "../src/objects/object_box_rotated.hpp"
#include "../src/objects/object_bump.hpp"
#include "../src/objects/object_curve.hpp"
#include "../src/objects/object_gate.hpp"
#include "../src/objects/object_ground.hpp"
#include "../src/objects/object_ground_rotated.hpp"
#include "../src/objects/object_ramp.hpp"
#include "../src/objects/object_silo.hpp"

#ifdef USE_VDB
#include "../src/objects/object_vdb.hpp"
#include "../src/sampling/sampling_particles_vdb.hpp"
#include "../src/sampling/regular_sampling_particles_vdb.hpp"
#endif

namespace py = pybind11;

// (Np, dim) row-major matrix used for bulk particle position/velocity access.
using ParticleArray = Eigen::Matrix<T, Eigen::Dynamic, TV::RowsAtCompileTime, Eigen::RowMajor>;

static ParticleArray toArray(const std::vector<TV>& vecs) {
    ParticleArray out(vecs.size(), TV::RowsAtCompileTime);
    for (size_t i = 0; i < vecs.size(); ++i)
        out.row(i) = vecs[i];
    return out;
}

static void fromArray(std::vector<TV>& vecs, const ParticleArray& arr) {
    if ((size_t)arr.rows() != vecs.size())
        throw std::runtime_error("Row count of array does not match the number of particles (Np).");
    for (size_t i = 0; i < vecs.size(); ++i)
        vecs[i] = arr.row(i);
}

PYBIND11_MODULE(matter, m) {
    m.doc() = "Python bindings for the matter MPM solver's Simulation class";

#ifdef USE_VDB
    openvdb::initialize();
#endif

    // ---------------- Enums ----------------
    py::enum_<PlateType>(m, "PlateType")
        .value("top", PlateType::top)
        .value("bottom", PlateType::bottom)
        .value("left", PlateType::left)
        .value("right", PlateType::right)
        .value("front", PlateType::front)
        .value("back", PlateType::back);

    py::enum_<ElasticModel>(m, "ElasticModel")
        .value("Hencky", ElasticModel::Hencky)
        .value("NeoHookean", ElasticModel::NeoHookean);

    py::enum_<PlasticModel>(m, "PlasticModel")
        .value("NoPlasticity", PlasticModel::NoPlasticity)
        .value("VM", PlasticModel::VM)
        .value("DP", PlasticModel::DP)
        .value("DPSoft", PlasticModel::DPSoft)
        .value("MCC", PlasticModel::MCC)
        .value("VMVisc", PlasticModel::VMVisc)
        .value("DPVisc", PlasticModel::DPVisc)
        .value("MCCVisc", PlasticModel::MCCVisc)
        .value("DPMui", PlasticModel::DPMui)
        .value("MCCMui", PlasticModel::MCCMui);

    py::enum_<HardeningLaw>(m, "HardeningLaw")
        .value("NoHard", HardeningLaw::NoHard)
        .value("ExpoExpl", HardeningLaw::ExpoExpl)
        .value("ExpoImpl", HardeningLaw::ExpoImpl)
        .value("SinhExpl", HardeningLaw::SinhExpl)
        .value("SinhImpl", HardeningLaw::SinhImpl);

    py::enum_<BC>(m, "BC")
        .value("NoSlip", BC::NoSlip)
        .value("SlipStick", BC::SlipStick)
        .value("SlipFree", BC::SlipFree);

    // ---------------- Particles ----------------
    py::class_<Particles>(m, "Particles")
        .def(py::init<unsigned int>(), py::arg("Np") = 1)
        .def_readwrite("x", &Particles::x)
        .def_readwrite("v", &Particles::v)
        .def_readwrite("pic", &Particles::pic)
        .def_readwrite("flip", &Particles::flip)
        .def_readwrite("eps_pl_dev", &Particles::eps_pl_dev)
        .def_readwrite("eps_pl_vol", &Particles::eps_pl_vol)
        .def_readwrite("eps_pl_vol_pradhana", &Particles::eps_pl_vol_pradhana)
        .def_readwrite("delta_gamma", &Particles::delta_gamma)
        .def_readwrite("viscosity", &Particles::viscosity)
        .def_readwrite("muI", &Particles::muI)
        .def_readwrite("F", &Particles::F)
        .def_readwrite("Bmat", &Particles::Bmat)
        .def_readwrite("tau", &Particles::tau)
        .def_readwrite("Ed", &Particles::Ed)
        .def_readwrite("minx_id", &Particles::minx_id)
        .def_readwrite("miny_id", &Particles::miny_id)
        .def_readwrite("minz_id", &Particles::minz_id)
        .def_readwrite("maxx_id", &Particles::maxx_id)
        .def_readwrite("maxy_id", &Particles::maxy_id)
        .def_readwrite("maxz_id", &Particles::maxz_id)
        .def("get_positions", [](const Particles& p) { return toArray(p.x); },
             "Returns particle positions as an (Np, dim) numpy array (a copy).")
        .def("set_positions", [](Particles& p, const ParticleArray& arr) { fromArray(p.x, arr); },
             "Sets particle positions from an (Np, dim) numpy array.")
        .def("get_velocities", [](const Particles& p) { return toArray(p.v); },
             "Returns particle velocities as an (Np, dim) numpy array (a copy).")
        .def("set_velocities", [](Particles& p, const ParticleArray& arr) { fromArray(p.v, arr); },
             "Sets particle velocities from an (Np, dim) numpy array.");

    // ---------------- Objects (ObjectGeneral hierarchy) ----------------
    py::class_<ObjectGeneral>(m, "ObjectGeneral")
        .def_readwrite("bc", &ObjectGeneral::bc)
        .def_readwrite("friction", &ObjectGeneral::friction)
        .def_readwrite("name", &ObjectGeneral::name)
        .def_readwrite("force_calc", &ObjectGeneral::force_calc)
        .def_readwrite("force", &ObjectGeneral::force);

    py::class_<ObjectBox, ObjectGeneral>(m, "ObjectBox")
        .def(py::init<BC, T, std::string, bool, TV, TV, TV>(),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0.0,
             py::arg("name") = "box", py::arg("force_calc") = false,
             py::arg("L") = TV::Ones(), py::arg("c") = TV::Zero(), py::arg("v") = TV::Zero())
        .def_readwrite("L", &ObjectBox::L)
        .def_readwrite("c", &ObjectBox::c)
        .def_readwrite("c0", &ObjectBox::c0)
        .def_readwrite("v", &ObjectBox::v);

    py::class_<ObjectBoxRotated, ObjectGeneral>(m, "ObjectBoxRotated")
        .def(py::init<BC, T, std::string, bool, TV, TV, TV, T>(),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0.0,
             py::arg("name") = "box_rotated", py::arg("force_calc") = false,
             py::arg("L") = TV::Ones(), py::arg("c") = TV::Zero(), py::arg("v") = TV::Zero(),
             py::arg("theta") = 0)
        .def_readwrite("L", &ObjectBoxRotated::L)
        .def_readwrite("c", &ObjectBoxRotated::c)
        .def_readwrite("c0", &ObjectBoxRotated::c0)
        .def_readwrite("v", &ObjectBoxRotated::v)
        .def_readwrite("ct", &ObjectBoxRotated::ct)
        .def_readwrite("st", &ObjectBoxRotated::st);

    py::class_<ObjectBump, ObjectGeneral>(m, "ObjectBump")
        .def(py::init<BC, T, std::string, bool>(),
             py::arg("bc"), py::arg("friction"),
             py::arg("name") = "bump", py::arg("force_calc") = false);

    py::class_<ObjectCurve, ObjectGeneral>(m, "ObjectCurve")
        .def(py::init<BC, T, std::string, bool>(),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0.0,
             py::arg("name") = "curve", py::arg("force_calc") = false);

    py::class_<ObjectGate, ObjectGeneral>(m, "ObjectGate")
        .def(py::init<BC, T, std::string, bool, T>(),
             py::arg("bc"), py::arg("friction"),
             py::arg("name") = "gate", py::arg("force_calc") = false, py::arg("height") = 0.016)
        .def_readwrite("height", &ObjectGate::height);

    py::class_<ObjectGround, ObjectGeneral>(m, "ObjectGround")
        .def(py::init<BC, T, std::string, bool, T>(),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0.0,
             py::arg("name") = "ground", py::arg("force_calc") = false, py::arg("y_ground") = 0)
        .def_readwrite("y_ground", &ObjectGround::y_ground);

    py::class_<ObjectGroundRotated, ObjectGeneral>(m, "ObjectGroundRotated")
        .def(py::init<BC, T, std::string, bool, T>(),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0.0,
             py::arg("name") = "", py::arg("force_calc") = false, py::arg("theta") = 0)
        .def_readwrite("ct", &ObjectGroundRotated::ct)
        .def_readwrite("st", &ObjectGroundRotated::st);

    py::class_<ObjectRamp, ObjectGeneral>(m, "ObjectRamp")
        .def(py::init<BC, T, std::string, bool>(),
             py::arg("bc"), py::arg("friction"),
             py::arg("name") = "ramp", py::arg("force_calc") = false);

#ifdef THREEDIM
    py::class_<ObjectSilo, ObjectGeneral>(m, "ObjectSilo")
        .def(py::init<BC, T, std::string, bool, T>(),
             py::arg("bc"), py::arg("friction"),
             py::arg("name") = "silo", py::arg("force_calc") = false, py::arg("cut") = -1)
        .def_readwrite("cut", &ObjectSilo::cut);
#endif

#ifdef USE_VDB
    py::class_<ObjectVdb, ObjectGeneral>(m, "ObjectVdb")
        .def(py::init<std::string, BC, T, std::string, bool>(),
             py::arg("filename"), py::arg("bc") = BC::NoSlip, py::arg("friction") = 0.0,
             py::arg("name") = "vdb", py::arg("force_calc") = false)
        .def("bounds", [](const ObjectVdb& self) {
            TV min_bbox = TV::Zero(), max_bbox = TV::Zero();
            self.bounds(min_bbox, max_bbox);
            return py::make_tuple(min_bbox, max_bbox);
        }, "Returns (min_bbox, max_bbox) of the levelset's active voxel bounding box.");
#endif

    // ---------------- ObjectPlate (standalone, not derived from ObjectGeneral) ----------------
#ifdef THREEDIM
    py::class_<ObjectPlate>(m, "ObjectPlate")
        .def(py::init<T, PlateType, BC, T, T, T, T, T, T, T, T, std::string, bool>(),
             py::arg("pos_object"), py::arg("plate_type"),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0,
             py::arg("pos_lower") = -1e15, py::arg("pos_upper") = 1e15,
             py::arg("vx_object") = 0.0, py::arg("vy_object") = 0.0, py::arg("vz_object") = 0.0,
             py::arg("vmin_factor") = 1.0, py::arg("load_factor") = 0.0,
             py::arg("name") = "plate", py::arg("force_calc") = false)
        .def_readwrite("vz_object", &ObjectPlate::vz_object)
        .def_readwrite("vz_object_original", &ObjectPlate::vz_object_original)
#else
    py::class_<ObjectPlate>(m, "ObjectPlate")
        .def(py::init<T, PlateType, BC, T, T, T, T, T, T, T, std::string, bool>(),
             py::arg("pos_object"), py::arg("plate_type"),
             py::arg("bc") = BC::NoSlip, py::arg("friction") = 0,
             py::arg("pos_lower") = -1e15, py::arg("pos_upper") = 1e15,
             py::arg("vx_object") = 0.0, py::arg("vy_object") = 0.0,
             py::arg("vmin_factor") = 1.0, py::arg("load_factor") = 0.0,
             py::arg("name") = "plate", py::arg("force_calc") = false)
#endif
        .def_readwrite("pos_object", &ObjectPlate::pos_object)
        .def_readwrite("plate_type", &ObjectPlate::plate_type)
        .def_readwrite("bc", &ObjectPlate::bc)
        .def_readwrite("friction", &ObjectPlate::friction)
        .def_readwrite("pos_upper", &ObjectPlate::pos_upper)
        .def_readwrite("pos_lower", &ObjectPlate::pos_lower)
        .def_readwrite("vx_object", &ObjectPlate::vx_object)
        .def_readwrite("vy_object", &ObjectPlate::vy_object)
        .def_readwrite("vx_object_original", &ObjectPlate::vx_object_original)
        .def_readwrite("vy_object_original", &ObjectPlate::vy_object_original)
        .def_readwrite("vmin_factor", &ObjectPlate::vmin_factor)
        .def_readwrite("load_factor", &ObjectPlate::load_factor)
        .def_readwrite("name", &ObjectPlate::name)
        .def_readwrite("force_calc", &ObjectPlate::force_calc)
        .def_readwrite("force", &ObjectPlate::force);

    // ---------------- Simulation ----------------
    py::class_<Simulation>(m, "Simulation")
        .def(py::init<>())
        .def_readonly("dim", &Simulation::dim)
        .def_readwrite("exit", &Simulation::exit)
        .def_readwrite("n_threads", &Simulation::n_threads)
        .def_readwrite("end_frame", &Simulation::end_frame)
        .def_readwrite("is_initialized", &Simulation::is_initialized)
        .def_readwrite("save_sim", &Simulation::save_sim)
        .def_readwrite("reduce_verbose", &Simulation::reduce_verbose)
        .def_readwrite("pbc", &Simulation::pbc)
        .def_readwrite("change_particle_positions", &Simulation::change_particle_positions)
        .def_readwrite("gravity_special", &Simulation::gravity_special)
        .def_readwrite("save_grid", &Simulation::save_grid)
        .def_readwrite("save_avg", &Simulation::save_avg)
        .def_readwrite("use_mibf", &Simulation::use_mibf)
        .def_readwrite("use_musl", &Simulation::use_musl)
        .def_readwrite("use_sparse", &Simulation::use_sparse)
        .def_readwrite("calculate_energy", &Simulation::calculate_energy)
        .def_readwrite("gravity", &Simulation::gravity)
        .def_readwrite("max_dt", &Simulation::max_dt)
        .def_readwrite("min_dt", &Simulation::min_dt)
        .def_readwrite("fps", &Simulation::fps)
        .def_readwrite("cfl", &Simulation::cfl)
        .def_readwrite("cfl_elastic", &Simulation::cfl_elastic)
        .def_readwrite("flip_ratio", &Simulation::flip_ratio)
        .def_readwrite("rho", &Simulation::rho)
        .def_readwrite("gravity_time", &Simulation::gravity_time)
        .def_readwrite("Lx", &Simulation::Lx)
        .def_readwrite("Ly", &Simulation::Ly)
#ifdef THREEDIM
        .def_readwrite("Lz", &Simulation::Lz)
#endif
        .def_readwrite("particles", &Simulation::particles)
        .def_readwrite("Np", &Simulation::Np)
        .def_readwrite("num_add_pbc_particles", &Simulation::num_add_pbc_particles)
        .def_readwrite("particle_mass", &Simulation::particle_mass)
        .def_readwrite("particle_volume", &Simulation::particle_volume)
        .def_readwrite("dx", &Simulation::dx)
        .def_readwrite("elastic_model", &Simulation::elastic_model)
        .def_readwrite("plastic_model", &Simulation::plastic_model)
        .def_readwrite("hardening_law", &Simulation::hardening_law)
        .def_readwrite("use_pradhana", &Simulation::use_pradhana)
        .def_readwrite("E", &Simulation::E)
        .def_readwrite("nu", &Simulation::nu)
        .def_readwrite("stress_tolerance", &Simulation::stress_tolerance)
        .def_readwrite("xi", &Simulation::xi)
        .def_readwrite("q_max", &Simulation::q_max)
        .def_readwrite("q_min", &Simulation::q_min)
        .def_readwrite("p_min", &Simulation::p_min)
        .def_readwrite("M", &Simulation::M)
        .def_readwrite("q_cohesion", &Simulation::q_cohesion)
        .def_readwrite("visc_exponent", &Simulation::visc_exponent)
        .def_readwrite("visc_time", &Simulation::visc_time)
        .def_readwrite("use_duvaut_lions_formulation", &Simulation::use_duvaut_lions_formulation)
        .def_readwrite("beta", &Simulation::beta)
        .def_readwrite("p0", &Simulation::p0)
        .def_readwrite("rho_s", &Simulation::rho_s)
        .def_readwrite("grain_diameter", &Simulation::grain_diameter)
        .def_readwrite("I_ref", &Simulation::I_ref)
        .def_readwrite("mu_1", &Simulation::mu_1)
        .def_readwrite("mu_2", &Simulation::mu_2)
        .def_readwrite("q_prefac", &Simulation::q_prefac)
        .def_readwrite("sampling_start_idx", &Simulation::sampling_start_idx)
        // pybind11 cannot steal a Python object into a std::unique_ptr<Base> function
        // parameter, so add_object is overloaded once per concrete derived type: each
        // overload copy-constructs a fresh C++ object (these classes hold only plain
        // value members, so the implicit copy constructor is safe) and wraps it in a
        // new unique_ptr owned by the simulation. Dispatch to the right overload is
        // based on the Python object's concrete registered type.
        .def("add_object", [](Simulation& self, const ObjectBox& obj) {
            self.objects.push_back(std::make_unique<ObjectBox>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectBoxRotated& obj) {
            self.objects.push_back(std::make_unique<ObjectBoxRotated>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectBump& obj) {
            self.objects.push_back(std::make_unique<ObjectBump>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectCurve& obj) {
            self.objects.push_back(std::make_unique<ObjectCurve>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectGate& obj) {
            self.objects.push_back(std::make_unique<ObjectGate>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectGround& obj) {
            self.objects.push_back(std::make_unique<ObjectGround>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectGroundRotated& obj) {
            self.objects.push_back(std::make_unique<ObjectGroundRotated>(obj));
        }, py::arg("object"))
        .def("add_object", [](Simulation& self, const ObjectRamp& obj) {
            self.objects.push_back(std::make_unique<ObjectRamp>(obj));
        }, py::arg("object"))
#ifdef THREEDIM
        .def("add_object", [](Simulation& self, const ObjectSilo& obj) {
            self.objects.push_back(std::make_unique<ObjectSilo>(obj));
        }, py::arg("object"))
#endif
#ifdef USE_VDB
        .def("add_object", [](Simulation& self, const ObjectVdb& obj) {
            self.objects.push_back(std::make_unique<ObjectVdb>(obj));
        }, py::arg("object"))
#endif
        .def("add_plate", [](Simulation& self, const ObjectPlate& obj) {
            self.plates.push_back(std::make_unique<ObjectPlate>(obj));
        }, py::arg("plate"), "Adds an ObjectPlate boundary (a copy of it, owned by the simulation).")
        .def_property_readonly("objects", [](Simulation& self) {
            std::vector<ObjectGeneral*> out;
            out.reserve(self.objects.size());
            for (auto& o : self.objects)
                out.push_back(o.get());
            return out;
        }, py::return_value_policy::reference_internal)
        .def_property_readonly("plates", [](Simulation& self) {
            std::vector<ObjectPlate*> out;
            out.reserve(self.plates.size());
            for (auto& o : self.plates)
                out.push_back(o.get());
            return out;
        }, py::return_value_policy::reference_internal)
        .def("initialize", &Simulation::initialize,
             py::arg("save") = true, py::arg("dir") = "output/", py::arg("name") = "dummy")
        .def("simulate", [](Simulation& self) {
            self.interrupt_check = []() { return PyErr_CheckSignals() != 0; };
            self.simulate();
            self.interrupt_check = nullptr;
            if (PyErr_Occurred())
                throw py::error_already_set();
            if (self.exit != 0)
                throw std::runtime_error("Simulation exited with an error (see console output above).");
        });

    // ---------------- Particle sampling ----------------
#ifdef THREEDIM
    m.def("sample_particles", [](Simulation& sim, T k_radius, T ppc, unsigned int crop_to_shape,
                                  std::uint32_t attempts, std::uint32_t seed) {
        sampleParticles(sim, k_radius, ppc, crop_to_shape, attempts, seed);
    }, py::arg("sim"), py::arg("k_radius"), py::arg("ppc") = 8, py::arg("crop_to_shape") = 0,
       py::arg("attempts") = 30, py::arg("seed") = 42);

    m.def("sample_particles_multi", [](Simulation& sim, std::vector<TV> origins, std::vector<TV> sizes,
                                        T k_radius, T ppc, std::uint32_t attempts, std::uint32_t seed) {
        sampleParticles(sim, origins, sizes, k_radius, ppc, attempts, seed);
    }, py::arg("sim"), py::arg("origins"), py::arg("sizes"), py::arg("k_radius"),
       py::arg("ppc") = 8, py::arg("attempts") = 30, py::arg("seed") = 42);
#else
    m.def("sample_particles", [](Simulation& sim, T k_radius, T ppc, unsigned int crop_to_shape,
                                  std::uint32_t attempts, std::uint32_t seed) {
        sampleParticles(sim, k_radius, ppc, crop_to_shape, attempts, seed);
    }, py::arg("sim"), py::arg("k_radius"), py::arg("ppc") = 6, py::arg("crop_to_shape") = 0,
       py::arg("attempts") = 200, py::arg("seed") = 42);

    m.def("sample_particles_multi", [](Simulation& sim, std::vector<TV> origins, std::vector<TV> sizes,
                                        T k_radius, T ppc, std::uint32_t attempts, std::uint32_t seed) {
        sampleParticles(sim, origins, sizes, k_radius, ppc, attempts, seed);
    }, py::arg("sim"), py::arg("origins"), py::arg("sizes"), py::arg("k_radius"),
       py::arg("ppc") = 6, py::arg("attempts") = 200, py::arg("seed") = 42);
#endif

    // regular_sample_particles/_multi place particles on a fixed grid instead of via
    // Poisson-disk sampling. Unlike sample_particles(...), sim.dx is an INPUT here
    // (used to derive particle spacing), not an output - set it before calling.
#ifdef THREEDIM
    m.def("regular_sample_particles", [](Simulation& sim, T ppc, unsigned int crop_to_shape) {
        regularSampleParticles(sim, ppc, crop_to_shape);
    }, py::arg("sim"), py::arg("ppc") = 8, py::arg("crop_to_shape") = 0);

    m.def("regular_sample_particles_multi", [](Simulation& sim, std::vector<TV> origins, std::vector<TV> sizes, T ppc) {
        regularSampleParticles(sim, origins, sizes, ppc);
    }, py::arg("sim"), py::arg("origins"), py::arg("sizes"), py::arg("ppc") = 8);
#else
    m.def("regular_sample_particles", [](Simulation& sim, T ppc, unsigned int crop_to_shape) {
        regularSampleParticles(sim, ppc, crop_to_shape);
    }, py::arg("sim"), py::arg("ppc") = 4, py::arg("crop_to_shape") = 0);

    m.def("regular_sample_particles_multi", [](Simulation& sim, std::vector<TV> origins, std::vector<TV> sizes, T ppc) {
        regularSampleParticles(sim, origins, sizes, ppc);
    }, py::arg("sim"), py::arg("origins"), py::arg("sizes"), py::arg("ppc") = 4);
#endif

#ifdef USE_VDB
#ifdef THREEDIM
    m.def("sample_particles_from_vdb", [](Simulation& sim, ObjectVdb& obj, T k_radius, T ppc) {
        sampleParticlesFromVdb(sim, obj, k_radius, ppc);
    }, py::arg("sim"), py::arg("obj"), py::arg("k_radius"), py::arg("ppc") = 8);

    m.def("sample_particles_from_vdb_multi", [](Simulation& sim, std::vector<ObjectVdb> objects, T k_radius, T ppc) {
        sampleParticlesFromVdb(sim, objects, k_radius, ppc);
    }, py::arg("sim"), py::arg("objects"), py::arg("k_radius"), py::arg("ppc") = 8);
#else
    m.def("sample_particles_from_vdb", [](Simulation& sim, ObjectVdb& obj, T k_radius, T ppc) {
        sampleParticlesFromVdb(sim, obj, k_radius, ppc);
    }, py::arg("sim"), py::arg("obj"), py::arg("k_radius"), py::arg("ppc") = 6);

    m.def("sample_particles_from_vdb_multi", [](Simulation& sim, std::vector<ObjectVdb> objects, T k_radius, T ppc) {
        sampleParticlesFromVdb(sim, objects, k_radius, ppc);
    }, py::arg("sim"), py::arg("objects"), py::arg("k_radius"), py::arg("ppc") = 6);
#endif

    // Same sim.dx-is-an-input caveat as regular_sample_particles above.
#ifdef THREEDIM
    m.def("regular_sample_particles_from_vdb", [](Simulation& sim, ObjectVdb& obj, T ppc) {
        regularSampleParticlesFromVdb(sim, obj, ppc);
    }, py::arg("sim"), py::arg("obj"), py::arg("ppc") = 8);

    m.def("regular_sample_particles_from_vdb_multi", [](Simulation& sim, std::vector<ObjectVdb> objects, T ppc) {
        regularSampleParticlesFromVdb(sim, objects, ppc);
    }, py::arg("sim"), py::arg("objects"), py::arg("ppc") = 8);
#else
    m.def("regular_sample_particles_from_vdb", [](Simulation& sim, ObjectVdb& obj, T ppc) {
        regularSampleParticlesFromVdb(sim, obj, ppc);
    }, py::arg("sim"), py::arg("obj"), py::arg("ppc") = 4);

    m.def("regular_sample_particles_from_vdb_multi", [](Simulation& sim, std::vector<ObjectVdb> objects, T ppc) {
        regularSampleParticlesFromVdb(sim, objects, ppc);
    }, py::arg("sim"), py::arg("objects"), py::arg("ppc") = 4);
#endif
#endif
}
