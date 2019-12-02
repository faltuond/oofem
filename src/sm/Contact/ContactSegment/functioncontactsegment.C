#include "functioncontactsegment.h"
#include "classfactory.h"



namespace oofem {

  REGISTER_ContactSegment(FunctionContactSegment);

    IRResultType FunctionContactSegment::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        
        /*int funcnum;
        IR_GIVE_FIELD(ir, funcnum, _IFT_FunctionContactSegment_function);*/

        //todo store function somehow

        //read circle parameters
        IR_GIVE_FIELD(ir, centerPoint, _IFT_FunctionContactSegment_centerpoint);
        IR_GIVE_FIELD(ir, radius, _IFT_FunctionContactSegment_radius);

        //todo check number of point coords matches domain

        return ContactSegment::initializeFrom(ir);
    }

    void FunctionContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tStep)
    {
        FloatArray nodeCoords;
        node->giveUpdatedCoordinates(nodeCoords, tStep);

        computeDistanceVector(answer, nodeCoords);

        //normalize
        double norm = answer.computeNorm();
        if ( norm > 1.0e-8 ) answer.times(1. / norm);
    }

    void FunctionContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        //returns just [[-1, 0], so that localization happens on node only
        //              [ 0,-1]]
        answer.resize(2, 2);
        answer.beUnitMatrix();
        answer.times(-1.); //?? seems reasonable to maintain compatibility with other segments
    }

    double FunctionContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        //get current nodal coords
        FloatArray nodeCoords, nodeCoordsInit, normal, normalInit;
        node->giveUpdatedCoordinates(nodeCoords, tStep);
        nodeCoordsInit = node->giveNodeCoordinates();

        computeDistanceVector(normal, nodeCoords);
        computeDistanceVector(normalInit, nodeCoordsInit);

        double cos = normal.dotProduct(normalInit);
        bool penetrated = cos <= 0;

        double answer = normal.computeNorm();
        if ( penetrated ) answer *= -1;

        return answer;
    }

    void FunctionContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
        s_loc.resize(0);
        //represents a function, does not have any dofs, so returns nothing
    }

    void FunctionContactSegment::computeDistanceVector(FloatArray & answer, const FloatArray & nodeCoords)
    {
        answer.beDifferenceOf(centerPoint, nodeCoords);
        double centerDistance = answer.computeNorm();
        answer.times((centerDistance - radius) / centerDistance);
    }

}
