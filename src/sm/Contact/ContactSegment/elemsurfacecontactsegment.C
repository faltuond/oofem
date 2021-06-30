#include "elemsurfacecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(Linear3dElementSurfaceContactSegment);



    void Linear3dElementSurfaceContactSegment::initializeFrom(InputRecord & ir)
    {
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

    void Linear3dElementSurfaceContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s) const
    {
    }

    void Linear3dElementSurfaceContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
    }

    bool Linear3dElementSurfaceContactSegment::computeContactPoint(FloatArray & ksi, Node * node, StructuralElement * element, int elemedge, TimeStep * tStep)
    {
        return false;
    }

}//end namespace OOFEM