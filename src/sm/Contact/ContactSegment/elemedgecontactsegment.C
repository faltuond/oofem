#include "elemedgecontactsegment.h"

namespace oofem {

    IRResultType ElementEdgeContactSegment::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        IR_GIVE_FIELD(ir, this->elemSet, _IFT_ElementEdgeContactSegment_elemSet);
        IR_GIVE_FIELD(ir, this->edgeSet, _IFT_ElementEdgeContactSegment_edgeSet);

        if ( elemSet.giveSize() != edgeSet.giveSize() ) {
            OOFEM_ERROR("The number of element indices and edge indices must correspond")
        }

        return ContactSegment::initializeFrom(ir);
    }

    void ElementEdgeContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep* tStep)
    {
        int nedges = edgeSet.giveSize();
        FloatArray normal, edgeVector, edgeNode1Coords, edgeNode2Coords, nodeCoords;
        double answerSize, normalSize;
        IntArray edgeNodes;

        node->giveUpdatedCoordinates(nodeCoords, tStep);

        //iterate over all edges, find the closest one
        for ( int pos = 0; pos < nedges; pos++ ) {

            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(elemSet(pos));
            element->giveBoundaryEdgeNodes(edgeNodes, edgeSet(pos));

            this->giveDomain()->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
            this->giveDomain()->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);

            edgeVector.beDifferenceOf(edgeNode2Coords, edgeNode1Coords);

            bool inbetween = computeDistanceVector(normal, nodeCoords, edgeNode1Coords, edgeNode2Coords);
            //no need to care here whether distance is negative or not
           
            double normalSize = sqrt(normal.dotProduct(normal));
            if (inbetween && (normalSize < answerSize || pos == 0 )) {
                answer = normal;
                answerSize = normalSize;
            }
        }
        //normalize
        answer.times(1. / answerSize);
    }

    void ElementEdgeContactSegment::computeExtendedNMatrix(FloatMatrix & answer, const Node * node)
    {
        
    }

    double ElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        int nedges = edgeSet.giveSize();
        FloatArray normal, edgeNode1Coords, edgeNode2Coords, nodeCoords;
        FloatArray normalInit, edgeNode1CoordsInit, edgeNode2CoordsInit, nodeCoordsInit;
        IntArray edgeNodes;
        double answerSize, normalSize;
        
        node->giveUpdatedCoordinates(nodeCoords, tStep);
        nodeCoordsInit = node->giveNodeCoordinates();

        //iterate over all edges, find the closest one
        for ( int pos = 0; pos < nedges; pos++ ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(elemSet(pos));
            element->giveBoundaryEdgeNodes(edgeNodes, edgeSet(pos));

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

            double normalSize = sqrt(normal.dotProduct(normal));
            if ( normalSize < answerSize ) {
                answerSize = normalSize;

                //todo work on this
            }
        }
        return answerSize;
    }

    void ElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {

    }

    bool ElementEdgeContactSegment::computeDistanceVector(FloatArray & answer, const FloatArray & externalPoint, const FloatArray & linePoint1, const FloatArray & linePoint2)
    {
        //this method has serious issues... for example when k = 0

        //convert edgevector to a line equation
        double k = (linePoint2(1) - linePoint1(1)) / (linePoint2(0) - linePoint1(0)); //slope
        double c1 = (linePoint1(1) - k * linePoint1(0)); //constant from y1 = k*x1 + c1

        //now we have the line equation in the form k*x + (-1)*y + c = 0
        //lets have a normal line in the form -1/k*x + (-1)*y + c2 = 0 and adjust c2 to make it pass through extern point
        double c2 = (1 / k)*externalPoint(0) + externalPoint(1);

        //find intersection point, for which y = k*x + c1 = (-1/k)*x + c2
        FloatArray isPoint(2);
        isPoint(0) = (c2 - c1) / (k + (1 / k));
        isPoint(1) = k * isPoint(0) + c1;

        answer.beDifferenceOf(isPoint, externalPoint);

        return abs(isPoint(0) - linePoint1(0)) < abs(linePoint2(0) - linePoint1(0)) && abs(isPoint(0) - linePoint2(0)) < abs(linePoint1(0) - linePoint2(0));
     }


}