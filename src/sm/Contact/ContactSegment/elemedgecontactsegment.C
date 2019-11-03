#include "elemedgecontactsegment.h"

namespace oofem {

    IRResultType ElementEdgeContactSegment::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        //IR_GIVE_FIELD(ir, this->elemSet, _IFT_ElementEdgeContactSegment_elemSet);
        int setnum;
        IR_GIVE_FIELD(ir, setnum, _IFT_ElementEdgeContactSegment_edgeSet);

        Set* set = this->giveDomain()->giveSet(setnum);
        if ( set == nullptr ) OOFEM_ERROR("Contact segment can not find set no. " + setnum);

        this->edges = set->giveBoundaryList();
        if ( edges.giveSize() <= 0 ) OOFEM_WARNING("Contact segment's edge list is empty");

        return ContactSegment::initializeFrom(ir);
    }

    void ElementEdgeContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep* tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //return zeros
            answer.resize(2);
            return;
        }

        StructuralElement* elem = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        int edgePos = closestEdge.at(2);

        FloatArray normal, edgeVector, edgeNode1Coords, edgeNode2Coords, nodeCoords;

        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes(edgeNodes, edgePos);

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        this->giveDomain()->giveNode(edgeNodes.at(1))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
        this->giveDomain()->giveNode(edgeNodes.at(2))->giveUpdatedCoordinates(edgeNode2Coords, tStep);

        bool inbetween = computeDistanceVector(normal, nodeCoords, edgeNode1Coords, edgeNode2Coords);
        //no need to care here whether distance is negative or not

        if ( inbetween ) {
            answer = normal;
            //normalize
            answer.times(1. / answer.computeNorm());
        }
        else {
            //todo maybe force reestablishing closest edge? instead of automatically giving up
            answer.resize(0);
        }
    }

    void ElementEdgeContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //return zeros
            answer.resize(6, 2);
            return;
        }

        StructuralElement* elem = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        int edgePos = closestEdge.at(2);
        
        FloatMatrix N;
        FloatArray contactPointCoords, dummyNormal, nodeCoords, edgeNode1Coords, edgeNode2Coords;
        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes(edgeNodes, edgePos);

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        this->giveDomain()->giveNode(edgeNodes.at(1))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
        this->giveDomain()->giveNode(edgeNodes.at(2))->giveUpdatedCoordinates(edgeNode2Coords, tStep);

        bool inbetween = computeDistanceVector(dummyNormal, nodeCoords, edgeNode1Coords, edgeNode2Coords, &contactPointCoords);

        //all the previous just to compute the contact point...

        elem->computeEdgeNMatrix(N, edgePos, contactPointCoords);

        answer.resize(N.giveNumberOfRows() + 2, N.giveNumberOfColumns());
        answer.zero();

        FloatMatrix extension(2, 2);
        extension.beUnitMatrix();

        answer.setSubMatrix(N, 1, 1);
        answer.setSubMatrix(extension, 1, N.giveNumberOfColumns() + 1);
    }

    double ElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //OOFEM_WARNING("Penetration asked of ContactSegment despite no contact occuring, returned 0.0");
            return 0.0;
        }

        FloatArray normal, edgeNode1Coords, edgeNode2Coords, nodeCoords;
        FloatArray normalInit, edgeNode1CoordsInit, edgeNode2CoordsInit, nodeCoordsInit;
        IntArray edgeNodes;

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        nodeCoordsInit = node->giveNodeCoordinates();

        StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        element->giveBoundaryEdgeNodes(edgeNodes, closestEdge.at(2));

        this->giveDomain()->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
        this->giveDomain()->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);
        edgeNode1CoordsInit = this->giveDomain()->giveNode(edgeNodes(0))->giveNodeCoordinates();
        edgeNode2CoordsInit = this->giveDomain()->giveNode(edgeNodes(1))->giveNodeCoordinates();

        bool inbetween = computeDistanceVector(normal, nodeCoords, edgeNode1Coords, edgeNode2Coords);
        bool inbetweenInit = computeDistanceVector(normalInit, nodeCoordsInit, edgeNode1CoordsInit, edgeNode2CoordsInit);

        //test whether initial and current vector are on different sides of line
        //find whether their cosine is larger than zero
        double cos = normal.dotProduct(normalInit);
        bool penetrated = cos <= 0;

        double answer = normal.computeNorm();
        if ( penetrated ) answer *= -1;

        return answer;
    }

    void ElementEdgeContactSegment::giveLocationArray( IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s)
    {
        answer.resize(0);
        //iterate over all segments
        for ( int pos = 0; pos < edges.giveSize() / 2; pos++ ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(edges(pos * 2));
            int edgePos = pos * 2 + 1;

            IntArray boundaryNodes;
            element->giveBoundaryEdgeNodes(boundaryNodes, edges(edgePos));

            IntArray locarray;
            element->giveBoundaryLocationArray(locarray, boundaryNodes, c_s, &dofIdArray);

            //add locarray to answer
            if ( pos == 0 ) {
                //in first iteration, try to guess future size of answer
                answer.resizeWithValues(answer.giveSize() + locarray.giveSize(), locarray.giveSize() * (edges.giveSize()/2));
            }
            else {
                //in other iterations just resize normally
                answer.resizeWithValues(answer.giveSize() + locarray.giveSize());
            }
            for ( int i = answer.giveSize() - locarray.giveSize(); i < answer.giveSize(); i++ ) answer.at(i) = locarray.at(i);
        }
    }

    void ElementEdgeContactSegment::giveClosestEdge(IntArray & answer, Node * node, TimeStep * tStep)
    {
        int knownIndex = giveIndexOfKnownNode(node);
        if ( knownIndex != -1 && knownClosestEdges.at(knownIndex).giveSize() == 2 ) {
            answer = knownClosestEdges.at(knownIndex);
            return;
        }
        //if previous if failed, it means that we don't know the closest edge to this node yet
        answer.resize(2);

        FloatArray normal, edgeNode1Coords, edgeNode2Coords, nodeCoords;
        IntArray edgeNodes;
        double answerSize = -1.;

        node->giveUpdatedCoordinates(nodeCoords, tStep);

        //iterate over all edges to find the closest one
        for ( int pos = 0; pos < edges.giveSize() / 2; pos++ ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(edges(pos * 2));
            int edgePos = pos * 2 + 1;

            element->giveBoundaryEdgeNodes(edgeNodes, edges(edgePos));

            this->giveDomain()->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
            this->giveDomain()->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);

            bool inbetween = computeDistanceVector(normal, nodeCoords, edgeNode1Coords, edgeNode2Coords);

            if ( inbetween ) {
                double normalSize = normal.computeNorm();
                if ( answerSize == -1. || normalSize < answerSize ) {
                    //new minimum found, update answer
                    answerSize = normalSize;
                    answer(0) = edges(pos * 2);
                    answer(1) = edges(edgePos);
                }
            }
        }
        //if answer size is still -1 it means no edge is properly in contact (no inbetween projections)
        if ( answerSize == -1. ) answer.resize(0);
        else {
            knownNodes.push_back(node);
            IntArray answerToStore(answer); //make copy to store to prevent it going out of scope
            knownClosestEdges.push_back(answerToStore);
        }
    }

    bool ElementEdgeContactSegment::computeDistanceVector(FloatArray & answer, const FloatArray & externalPoint, const FloatArray & linePoint1, const FloatArray & linePoint2, /*optional*/ FloatArray * contactPointCoords)
    {

        //convert edgevector to an equation
        FloatArray lineVector;
        lineVector.beDifferenceOf(linePoint2, linePoint1);

        //if we consider the normal line in the form ax + by + c = 0, an and bn are components of lineVector
        //(which is normal to normal line)
        double an = lineVector.at(1);
        double bn = lineVector.at(2);
        //c has to be found using the external point coordinates
        double cn = -an * externalPoint.at(1) - bn * externalPoint.at(2);

        //now create an vector that is parallel to the normal (i.e. perpendicular to the edge line)
        FloatArray normalVector(2);
        normalVector.at(1) = -lineVector.at(2);
        normalVector.at(2) = lineVector.at(1);

        //thus we have equation of the edge line
        double ae = lineVector.at(1);
        double be = lineVector.at(2);
        double ce = -ae * linePoint1.at(1) - be * linePoint1.at(2);

        //now we have to find the intersecting point
        double divider = ae * bn - an * be;//should only be zero if lines are parallel
        FloatArray contactPoint(2);
        contactPoint.at(1) = (ce*bn - cn*be) / divider;
        contactPoint.at(2) = (ae*cn - an*ce) / divider;

        answer.beDifferenceOf(contactPoint, externalPoint);

        if ( contactPointCoords != nullptr ) {
            contactPointCoords->resize(2);
            contactPointCoords->at(1) = contactPoint.at(1);
            contactPointCoords->at(2) = contactPoint.at(2);
        }

        double lineLength = abs(linePoint1(0) - linePoint2(0));
        //point lies inbetween line points if it is closer to both than length of line
        return abs(contactPoint(0) - linePoint1(0)) < lineLength && abs(contactPoint(0) - linePoint2(0)) < lineLength;
    }

void
ElementEdgeContactSegment :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
  // do the update
}


  
}


