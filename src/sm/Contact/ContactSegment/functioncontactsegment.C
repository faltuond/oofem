#include "functioncontactsegment.h"
#include "classfactory.h"



namespace oofem {

  //REGISTER_ContactSegment(FunctionContactSegment);

   
    void FunctionContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tStep)
    {
        FloatArray nodeCoords;
        node->giveUpdatedCoordinates(nodeCoords, tStep);

        computeDistanceVector(answer, nodeCoords);

        //normalize according to normalization mode specified
        if ( normmode != NM_Never ) {
            double norm = answer.computeNorm();
            if ( normmode == NM_Always || norm > 1.0e-8 ) answer.times(1. / norm);
        }
    }

    void FunctionContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        //returns just [[-1, 0], so that localization happens on node only
        //              [ 0,-1]]
        //adapted to size of node coordinates to make the class independent on dimension
        int ncoords = node->giveCoordinates()->giveSize();
        answer.resize(ncoords, ncoords);
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

    void FunctionContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
        s_loc.resize(0);
    }

}
