#include "elemsurfacecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(Linear3dElementSurfaceContactSegment);



    void Linear3dElementSurfaceContactSegment::initializeFrom(InputRecord & ir)
    {
        BoundaryContactSegment::initializeFrom(ir);
    }

    void Linear3dElementSurfaceContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tStep)
    {
        if (hasNonLinearGeometry(node, tStep)) {
            OOFEM_ERROR("Nonlinear 3D contact segments not yet implemented");
        }
        else {
            //determine normal from surface
            IntArray closestSurface;
            giveClosestBoundary(closestSurface, node, tStep);
            if (closestSurface.giveSize() != 2) {
                //no closest surface means no contact
                //return zeros
                answer.resize(3);
                return;
            }

            StructuralElement *elem = (StructuralElement *)this->giveDomain()->giveElement(closestSurface.at(1));
            int surfacePos = closestSurface.at(2);

            FloatArray cPointLocal, normal;

            bool inbetween = computeContactPoint(cPointLocal, node, elem, surfacePos, tStep);
            //no need to care here whether distance is negative or not

            //retrieve edge normal from element interpolation
            FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>(elem->giveInterpolation());
            if (interpolation == nullptr) {
                OOFEM_ERROR("Non-3D element encountered in Linear3dElementSurfaceContactSegment");
            }
            interpolation->surfaceEvalNormal(normal, surfacePos, cPointLocal, FEIElementGeometryWrapper(elem));

            answer = normal;
            answer.times(-1.);
        }
    }

    void Linear3dElementSurfaceContactSegment::computeTangent(FloatArray & answer, Node * node, TimeStep * tstep)
    {
        OOFEM_ERROR("Nonlinear 3D contact segments not yet implemented");
    }

    void Linear3dElementSurfaceContactSegment::computeSegmentNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestBoundary(closestEdge, node, tStep);
        if (closestEdge.giveSize() != 2) {
            //no closest edge means no contact
            //return zeros
            answer.resize(2, 4);
            return;
        }

        StructuralElement *elem = (StructuralElement *)this->giveDomain()->giveElement(closestEdge.at(1));
        int surfacePos = closestEdge.at(2);

        FloatMatrix N;
        FloatArray cPointLocal;

        bool inbetween = computeContactPoint(cPointLocal, node, elem, surfacePos, tStep);
        //all the previous just to compute the contact point...
        elem->computeEdgeNMatrix(N, surfacePos, cPointLocal);
    }

    void Linear3dElementSurfaceContactSegment::computeSegmentBMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        OOFEM_ERROR("Nonlinear 3D contact segments not yet implemented");
    }

    double Linear3dElementSurfaceContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        IntArray closestSurface;
        giveClosestBoundary(closestSurface, node, tStep);
        if (closestSurface.giveSize() != 2) {
            //no closest surface means no contact
            return 0.0;
        }

        FloatArray cPointLocal, projection, nodeCoords, surfaceNode1Coords;
        IntArray surfaceNodes;

        node->giveUpdatedCoordinates(nodeCoords, tStep);


        StructuralElement *element = (StructuralElement *)this->giveDomain()->giveElement(closestSurface.at(1));

        surfaceNodes = element->giveBoundarySurfaceNodes(closestSurface.at(2));
        element->giveNode(surfaceNodes(0))->giveUpdatedCoordinates(surfaceNode1Coords, tStep);

        bool inbetween = computeContactPoint(cPointLocal, node, element, closestSurface.at(2), tStep);
        if (inbetween == false) {
            return 0;
        }
        //projection.beDifferenceOf(nodeCoords, cPoint);
        projection.beDifferenceOf(nodeCoords, surfaceNode1Coords);

        ////test whether initial and current vector are on different sides of line
        FloatArray normal;
        this->computeNormal(normal, node, tStep);
        return projection.dotProduct(normal);
    }

    void Linear3dElementSurfaceContactSegment::computeMetricTensor(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        OOFEM_ERROR("Nonlinear 3D contact segments not yet implemented");
    }

    void Linear3dElementSurfaceContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s) const
    {
        if (lastBoundary.giveSize() == 2) {
            StructuralElement *element = (StructuralElement *)this->giveDomain()->giveElement(lastBoundary.at(1));
            int boundaryPos = lastBoundary.at(2);

            IntArray boundaryNodes = element->giveBoundaryEdgeNodes(boundaryPos);
            element->giveBoundaryLocationArray(answer, boundaryNodes, c_s, nullptr);
        }
        else {
            //if no segment was worked with, returns zeros
            //@todo What array size to return? How do we know how many nodes per boundary are there when we do not have any boundary? (2D case returns 4 here)
            answer.resize(0);
        }
    }
    
    //@todo This might probably be moved to parent class, as it is precisely the same as in 2D
    void Linear3dElementSurfaceContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s)
    {
        answer.resize(0);
        IntArray boundaryloc, boundaryNodes;
        //iterate over all edges and add their locarrays
        for (int pos = 0; pos < boundaries.giveSize() / 2; pos++) {
            StructuralElement *element = (StructuralElement *)this->giveDomain()->giveElement(boundaries(pos * 2));
            int boundaryPos = pos * 2 + 1;

            element->giveBoundaryLocationArray(boundaryloc, element->giveBoundaryEdgeNodes(boundaries(boundaryPos)), c_s, nullptr);

            answer.followedBy(boundaryloc);
        }
    }

    bool Linear3dElementSurfaceContactSegment::computeContactPoint(FloatArray & ksi, Node * node, StructuralElement * elem, int elemsurface, TimeStep * tStep)
    {
        //@todo Projection on surface, returning parametric coordinates ksi

        FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>(elem->giveInterpolation());
        if (interpolation == nullptr) {
            OOFEM_ERROR("Non-3D element encountered in Linear3dElementSurfaceContactSegment");
        }

        //This is a 3D linear element, therefore its boundary is a surface with at least 3 nodes, of which no 3 are in a straight line
        FloatArray extNodeCoords, surfNodeCoords1, surfNodeCoords2, surfNodeCoords3;
        IntArray surfaceNodes;

        surfaceNodes = elem->giveBoundarySurfaceNodes(elemsurface);
        node->giveUpdatedCoordinates(extNodeCoords, tStep);
        elem->giveNode(surfaceNodes(0))->giveUpdatedCoordinates(surfNodeCoords1, tStep);
        elem->giveNode(surfaceNodes(1))->giveUpdatedCoordinates(surfNodeCoords2, tStep);
        elem->giveNode(surfaceNodes(2))->giveUpdatedCoordinates(surfNodeCoords3, tStep);

        FloatArray tangent1, tangent2, normal;

        tangent1.beDifferenceOf(surfNodeCoords2, surfNodeCoords1);
        tangent2.beDifferenceOf(surfNodeCoords3, surfNodeCoords1);
        normal.beVectorProductOf(tangent1, tangent2);
        normal.normalize();

        FloatArray projection, normalProjection;
        projection.beDifferenceOf(extNodeCoords, surfNodeCoords1);
        double distance = normal.dotProduct(projection);
        normalProjection = normal;
        normalProjection.times(distance);

        FloatArray contactPoint;
        contactPoint.beDifferenceOf(extNodeCoords, normalProjection);

        int inbetween = interpolation->global2local(ksi, contactPoint, FEIElementDeformedGeometryWrapper(elem));
        if (inbetween == 0) {
            //outside of surface, return zeros and false
            ksi.resize(2);
            return false;
        }

        //We obtained the parametric coordinates for the whole element. We need them only for the surface
        //We are reasonably sure, though, that one of those is 1 or -1, since we are on the surface of the element
        //That one should be eliminated.

        //NO. What if we are exactly in the vertex and all of them are 1?




        return true;
    }

}//end namespace OOFEM