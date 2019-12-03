
#include "qfunctioncontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(QFunctionContactSegment)

    IRResultType QFunctionContactSegment::initializeFrom(InputRecord* ir) {
        IRResultType result;

        IR_GIVE_FIELD(ir, a, _IFT_QFunctionContactSegment_a);
        IR_GIVE_FIELD(ir, b, _IFT_QFunctionContactSegment_b);
        IR_GIVE_FIELD(ir, c, _IFT_QFunctionContactSegment_c);

        return FunctionContactSegment::initializeFrom(ir);

    }

}