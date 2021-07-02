#include "elemsurfacecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(Linear3dElementSurfaceContactSegment);



    void Linear3dElementSurfaceContactSegment::initializeFrom(InputRecord & ir)
    {
        BoundaryContactSegment::initializeFrom(ir);
    }

    void Linear3dElementSurfaceContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tstep)
    {

    }

    void Linear3dElementSurfaceContactSegment::computeTangent(FloatArray & answer, Node * node, TimeStep * tstep)
    {
    }

    void Linear3dElementSurfaceContactSegment::computeSegmentNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
    }

    void Linear3dElementSurfaceContactSegment::computeSegmentBMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
    }

    double Linear3dElementSurfaceContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        return 0.0;
    }

    void Linear3dElementSurfaceContactSegment::computeMetricTensor(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
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

    bool Linear3dElementSurfaceContactSegment::computeContactPoint(FloatArray & ksi, Node * node, StructuralElement * element, int elemedge, TimeStep * tStep)
    {
        //@todo Projection on surface, returning parametric coordinates ksi
        return false;
    }

}//end namespace OOFEM