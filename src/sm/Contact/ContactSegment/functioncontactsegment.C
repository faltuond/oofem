#include "functioncontactsegment.h"
#include "classfactory.h"



namespace oofem {

  //REGISTER_ContactSegment(FunctionContactSegment);

   
    void FunctionContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tStep)
    {
        FloatArray nodeCoords;
        node->giveUpdatedCoordinates(nodeCoords, tStep);

        FloatArray dummyContactPoint;
        computeContactPoint(dummyContactPoint, answer, nodeCoords); //should already be normalized

        //normalize according to normalization mode specified
        /*if ( normmode != NM_Never ) {
            double norm = answer.computeNorm();
            if ( normmode == NM_Always || norm > 1.0e-8 ) answer.times(1. / norm);
        }*/
    }

	void FunctionContactSegment::computeTangent( FloatArray &answer, Node *node, TimeStep *tstep )
    {
        FloatArray normal;
        computeNormal( normal, node, tstep );
        answer.resize( 2 ); //this is NOT a strictly 2D segment
        /*answer.at( 1 ) = normal.at( 2 );
        answer.at( 2 ) = -normal.at( 1 );*/

		OOFEM_ERROR( "Tangent vector for function segments not implemented" );
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

    void FunctionContactSegment::computeExtendedBMatrix( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
		//is only the nodal part,same as extN, which in this case is all zeros
        int ncoords = node->giveCoordinates()->giveSize();
        answer.resize( ncoords, ncoords );
    }

    bool FunctionContactSegment::hasNonLinearGeometry( Node *node, TimeStep *tStep )
    {
        return false; //placeholder??
    }

    void FunctionContactSegment::computeMetricTensor( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
        OOFEM_ERROR( "Not implemented for function segments" );
    }

    double FunctionContactSegment::computePenetration( Node *node, TimeStep *tStep )
    {
        //get current nodal coords
        FloatArray nodeCoords, nodeCoordsInit, contactPoint, contactPointInit, dummyNormal;
        node->giveUpdatedCoordinates(nodeCoords, tStep);
        nodeCoordsInit = node->giveNodeCoordinates();

        computeContactPoint(contactPoint, dummyNormal, nodeCoords);
        computeContactPoint(contactPointInit, dummyNormal, nodeCoordsInit);

        FloatArray projection, projectionInit;
        projection.beDifferenceOf(contactPoint, nodeCoords);
        projectionInit.beDifferenceOf(contactPointInit, nodeCoordsInit);

        double cos = projection.dotProduct(projectionInit);
        bool penetrated = cos <= 0;

        double answer = projection.computeNorm();
        if ( penetrated ) answer *= -1;

        return answer;
    }


    void FunctionContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s) const
    {
        s_loc.resize(0);
        //represents a function, does not have any dofs, so returns nothing
    }

    void FunctionContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
        s_loc.resize(0);
    }

}
