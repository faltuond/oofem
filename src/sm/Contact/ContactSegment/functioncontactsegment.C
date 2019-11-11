#include "functioncontactsegment.h"

namespace oofem {
       
    IRResultType FunctionContactSegment::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        
        int funcnum;
        IR_GIVE_FIELD(ir, funcnum, _IFT_FunctionContactSegment_function);

        //todo store function somehow

        return ContactSegment::initializeFrom(ir);
    }

    void FunctionContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tstep)
    {
        //todo
    }

    void FunctionContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        //returns just [[1,0], so that localization happens on node only
        //              [0,1]]
        answer.resize(2, 2);
        answer.beUnitMatrix();
    }

    double FunctionContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        //todo
        return 0.0;
    }

    void FunctionContactSegment::giveLocationArray(IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
        s_loc.resize(0);
        //represents a function, does not have any dofs
    }

    void FunctionContactSegment::updateYourself(TimeStep * tStep)
    {
        //may be unnecessary
    }

}