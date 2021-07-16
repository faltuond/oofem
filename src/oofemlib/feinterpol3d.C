/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "feinterpol3d.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
double FEInterpolation3d :: giveVolume(const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented in subclass.");
    return 0;
}

IntArray FEInterpolation3d :: boundaryEdgeGiveNodes(int boundary) const
{
    return this->computeLocalEdgeMapping(boundary);
}

void FEInterpolation3d :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

IntArray FEInterpolation3d :: boundaryGiveNodes(int boundary) const
{
    return this->computeLocalSurfaceMapping(boundary);
}

void FEInterpolation3d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->surfaceEvalNormal(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->surfaceGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->surfaceLocal2global(answer, boundary, lcoords, cellgeo);
}

IntArray FEInterpolation3d :: computeEdgeMapping(const IntArray &elemNodes, int iedge) const
{
    const auto &ln = this->computeLocalEdgeMapping(iedge);
    int size = ln.giveSize();
    IntArray edgeNodes(size);
    for ( int i = 1; i <= size; i++ ) {
        edgeNodes.at(i) = elemNodes.at( ln.at(i) );
    }
    return edgeNodes;
}

#define POINT_TOL 1.e-3

int FEInterpolation3d::surfaceGlobal2local(FloatArray & answer, int isurf, const FloatArray & gcoords, const FEICellGeometry & cellgeo)
{
    //What follows is a copy of FEInterpol2d::global2local(), just with the functions invoked changed to fit the case with 3D element's surface

    FloatArray res, delta, delta_3d, guess, lcoords_guess, lcoords_guess_3d;
    IntArray surface_coord_indices;
    FloatMatrix jac;
    double convergence_limit, error = 0.0;
    double const_surface_coord;

    // depending on which surface we are on, find the indices of the coords relevant to the surface
    // and also the value of the one which is irrelevant
    const_surface_coord = this->surfaceGiveLCoordIndices(surface_coord_indices, isurf);

    // find a suitable convergence limit
    convergence_limit = 1e-6 /* this->giveCharacteristicLength(cellgeo)*/; //@todo include char. length of surface? How to get it?

    // setup initial guess
    lcoords_guess.resize(2);
    lcoords_guess.zero();

    // apply Newton-Raphson to solve the problem
    for (int nite = 0; nite < 10; nite++) {
        // compute the residual
        this->surfaceLocal2global(guess, isurf, lcoords_guess, cellgeo);
        res = { gcoords(0) - guess(0), gcoords(1) - guess(1), gcoords(2) - guess(2) };

        //@todo create 3D guess based on what surface we are on
        lcoords_guess_3d.resize(3);
        lcoords_guess_3d.add(const_surface_coord);
        lcoords_guess_3d.at(surface_coord_indices.at(1)) = lcoords_guess.at(1);
        lcoords_guess_3d.at(surface_coord_indices.at(2)) = lcoords_guess.at(2);

        // check for convergence
        error = res.computeNorm();
        if (error < convergence_limit) {
            break;
        }

        // compute the corrections
        this->surfaceGiveJacobianMatrixAt(jac, isurf, lcoords_guess, cellgeo);
        jac.solveForRhs(res, delta_3d);

        //reduce 3d guess back to 2d
        delta.resize(2);
        delta.at(1) = delta_3d.at(surface_coord_indices.at(1));
        delta.at(2) = delta_3d.at(surface_coord_indices.at(2));

        // update guess
        lcoords_guess.add(delta);

    }
    if (error > convergence_limit) { // Imperfect, could give false negatives.
        OOFEM_WARNING("Failed convergence");
        answer = { 1. / 3., 1. / 3. };
        return false;
    }

    answer = { lcoords_guess(0), lcoords_guess(1) };

    // test if inside
    bool inside = true;
    for (int i = 1; i <= 2; i++) {
        if (answer.at(i) < (-1. - POINT_TOL)) {
            answer.at(i) = -1.;
            inside = false;
        }
        else if (answer.at(i) > (1. + POINT_TOL)) {
            answer.at(i) = 1.;
            inside = false;
        }
    }

    return inside;
}

void FEInterpolation3d::surfaceGiveJacobianMatrixAt(FloatMatrix & jacobianMatrix, int isurf, const FloatArray & lcoords, const FEICellGeometry & cellgeo)
{
    OOFEM_ERROR("Surface Jacobian matrix: Not implemented in a general way, needs to be overloaded if desired");
}

void FEInterpolation3d::surfaceEvalBaseVectorsAt(FloatArray & G1, FloatArray & G2, int isurf, const FloatArray & lcoords, const FEICellGeometry & cellgeo)
{
    OOFEM_ERROR("Not implemented");
}

double FEInterpolation3d::surfaceGiveLCoordIndices(IntArray &answer, const int isurf) const {
    OOFEM_ERROR("Not implemented");
    return 0.0;
}

IntArray FEInterpolation3d :: computeSurfaceMapping(const IntArray &elemNodes, int isurf) const
{
    const auto &ln = this->computeLocalSurfaceMapping(isurf);
    int size = ln.giveSize();
    IntArray surfNodes(size);
    for ( int i = 1; i <= size; i++ ) {
        surfNodes.at(i) = elemNodes.at( ln.at(i) );
    }
    return surfNodes;
}

std::unique_ptr<IntegrationRule> FEInterpolation3d :: giveBoundaryEdgeIntegrationRule(int order, int boundary)
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + this->order);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return iRule;
}

void FEInterpolation3d :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not implemented");
}

void FEInterpolation3d :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not implemented");
}

void FEInterpolation3d::surfaceEvaldNdxi(FloatMatrix & answer, int isurf, const FloatArray & lcoords, const FEICellGeometry & cellgeo)
{
    OOFEM_ERROR("Not implemented");
}

double FEInterpolation3d :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not implemented");
    return -1.0;
}

IntArray FEInterpolation3d::boundarySurfaceGiveNodes(int boundary) const
{
    return this->computeLocalSurfaceMapping(boundary);
}
  
} // end namespace oofem
