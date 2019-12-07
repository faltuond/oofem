#include "circlecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(CircleContactSegment);

    IRResultType CircleContactSegment::initializeFrom(InputRecord* ir) {
        //read circle parameters
        IRResultType result;
        IR_GIVE_FIELD(ir, centerPoint, _IFT_CircleContactSegment_centerpoint);
        IR_GIVE_FIELD(ir, radius, _IFT_CircleContactSegment_radius);

        return FunctionContactSegment::initializeFrom(ir);
    }

    void CircleContactSegment::computeDistanceVector(FloatArray & answer, const FloatArray & nodeCoords)
    {
        if ( nodeCoords.giveSize() != centerPoint.giveSize() ) {
            OOFEM_ERROR("Node coordinate dimension (%i) does not match circle/sphere center point dimension (%i)", nodeCoords.giveSize(), centerPoint.giveSize());
        }
        answer.beDifferenceOf(centerPoint, nodeCoords);
        double centerDistance = answer.computeNorm();
        answer.times((centerDistance - radius) / centerDistance);
    }
}//end namespace oofem