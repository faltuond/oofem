#include "circlecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(CircleContactSegment);

    IRResultType CircleContactSegment::initializeFrom(InputRecord* ir) {
        //read circle parameters
        IRResultType result;
        IR_GIVE_FIELD(ir, centerPoint, _IFT_CircleContactSegment_centerpoint);
        IR_GIVE_FIELD(ir, radius, _IFT_CircleContactSegment_radius);

        //todo check number of point coords matches domain

        return FunctionContactSegment::initializeFrom(ir);
    }

    void CircleContactSegment::computeDistanceVector(FloatArray & answer, const FloatArray & nodeCoords)
    {
        answer.beDifferenceOf(centerPoint, nodeCoords);
        double centerDistance = answer.computeNorm();
        answer.times((centerDistance - radius) / centerDistance);
    }
}//end namespace oofem