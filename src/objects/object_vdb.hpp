// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef VDB_OBJ_H
#define VDB_OBJ_H

#include "object_general.hpp"

#define IMATH_HALF_NO_LOOKUP_TABLE
#include <openvdb/openvdb.h>
#include <openvdb/io/File.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/Interpolation.h>

// Make sure to fill the interor when creating the levelset, otherwise normals
// and signed distances are not computed correctly

class ObjectVdb : public ObjectGeneral {
public:
    typedef typename openvdb::Grid<typename openvdb::tree::Tree4<float, 5, 4, 3>::Type> GridT;
    typedef typename GridT::TreeType TreeT;
    typedef typename openvdb::tools::ScalarToVectorConverter<GridT>::Type GradientGridT;
    typedef typename GradientGridT::TreeType GradientTreeT;

    typename GridT::Ptr grid;
    typename GradientGridT::Ptr grad_phi;

    ~ObjectVdb(){}

    ObjectVdb(std::string filename, BC bc_in = BC::NoSlip, T friction_in = 0.0, std::string name_in = "vdb", bool force_calc_in = false) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in) {

        openvdb::io::File file(filename);
        file.open();
        openvdb::GridPtrVecPtr my_grids = file.getGrids();
        file.close();
        int count = 0;
        for (openvdb::GridPtrVec::iterator iter = my_grids->begin(); iter != my_grids->end(); ++iter) {
            grid = openvdb::gridPtrCast<GridT>(*iter);
            count++;
        }

        openvdb::tools::Gradient<GridT> mg(*grid);
        grad_phi = mg.process();
    }

    bool inside(const TV& X_in) const override {
        int dim = X_in.size();
        Eigen::Matrix<T, 3, 1> X;
        X.setZero();
        for (int d = 0; d < dim; d++)
            X(d) = X_in(d);

        openvdb::tools::GridSampler<TreeT, openvdb::tools::BoxSampler> interpolator(grid->constTree(), grid->transform());
        openvdb::math::Vec3<T> P(X(0), X(1), X(2));
        float phi = interpolator.wsSample(P); // this is the signed distance

        return ((T)phi <= 0);
    }

    TV normal(const TV& X_in) const override {
        int dim = X_in.size();
        Eigen::Matrix<T, 3, 1> X;
        X.setZero();
        for (int d = 0; d < dim; d++)
            X(d) = X_in(d);

        openvdb::tools::GridSampler<GradientTreeT, openvdb::tools::BoxSampler> interpolator(grad_phi->constTree(), grad_phi->transform());
        openvdb::math::Vec3<T> P(X(0), X(1), X(2));
        auto grad_phi = interpolator.wsSample(P);
        TV result;
        for (int d = 0; d < dim; d++)
            result(d) = grad_phi(d);
        T norm = result.norm();
        if (norm != 0)
            return result / norm;
        else
            return TV::Zero();
    }

    void bounds(TV& min_bbox, TV& max_bbox) const {
        int dim = min_bbox.size();
        min_bbox.setZero();
        max_bbox.setZero();
        openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
        auto wmin = grid->indexToWorld(bbox.min());
        auto wmax = grid->indexToWorld(bbox.max());

        for (int d = 0; d < dim; d++) {
            min_bbox(d) = (T)wmin(d);
            max_bbox(d) = (T)wmax(d);
        }
    }


}; // End class ObjectVdb

#endif // VDB_OBJ_H
