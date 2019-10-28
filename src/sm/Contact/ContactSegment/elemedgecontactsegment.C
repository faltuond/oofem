#include "elemedgecontactsegment.h"

namespace oofem {

    IRResultType ElementEdgeContactSegment::initializeFrom(InputRecord * Air)
    {
        
    }

    void ElementEdgeContactSegment::computeNormal(FloatArray & answer, const Node * node)
    {
    }

    void ElementEdgeContactSegment::computeExtendedNMatrix(FloatMatrix & answer, const Node * node)
    {
    }

    double ElementEdgeContactSegment::computePenetration(const Node * node)
    {
        return 0.0;
    }

    void ElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
    }


}