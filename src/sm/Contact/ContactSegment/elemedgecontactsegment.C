#include "elemedgecontactsegment.h"
#include "classfactory.h"

namespace oofem {


    REGISTER_ContactSegment(ElementEdgeContactSegment);


    IRResultType ElementEdgeContactSegment::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        //IR_GIVE_FIELD(ir, this->elemSet, _IFT_ElementEdgeContactSegment_elemSet);
        IR_GIVE_FIELD(ir, setnum, _IFT_ElementEdgeContactSegment_edgeSet);

        updateMode = UM_EachStep;
        int updateModeInt = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, updateModeInt, _IFT_ElementEdgeContactSegment_pairUpdateMode);
        if ( result == IRRT_OK ) {
            updateMode = (UpdateMode)updateModeInt;
            if ( updateModeInt < 0 || updateModeInt > 2 ) OOFEM_ERROR("Contact segment pairing update mode can be only 0, 1 or 2");
        }

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

        FloatArray cPoint, cPointLocal, normal, edgeVector, edgeNode1Coords, edgeNode2Coords, nodeCoords;

        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes(edgeNodes, edgePos);

        node->giveUpdatedCoordinates(nodeCoords, tStep);

        elem->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
        elem->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);

        bool inbetween = computeContactPoint(cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords);
        //no need to care here whether distance is negative or not

        //retrieve edge normal from element interpolation
        FEInterpolation2d* interpolation = dynamic_cast<FEInterpolation2d*>(elem->giveInterpolation());
        if ( interpolation == nullptr ) {
            OOFEM_ERROR("Non-2D element encountered in ElementEdgeContactSegment");
        }

        interpolation->global2local(cPointLocal, cPoint, FEIElementGeometryWrapper(elem));
        interpolation->edgeEvalNormal(normal, edgePos, cPointLocal, FEIElementGeometryWrapper(elem));

        //if element geometry is nonlinear, normal should be converted to deformed configuration
        NLStructuralElement* nlelem = dynamic_cast<NLStructuralElement*>(elem);
        if ( nlelem && nlelem->giveGeometryMode() == 1 ) transformNormalToDeformedShape(normal, nlelem, cPointLocal, tStep);


        //contact normal is defined as going towards the edge
        normal.times(-1.);
        answer = normal;


        //if ( true || inbetween ) { //attempt to fix concave problem
        //    answer = normal;
        //    //normalize according to normalization mode specified
        //    if ( normmode != NM_Never ) {
        //        double norm = answer.computeNorm();
        //        if ( (normmode == NM_Always || norm > 1.0e-8) && norm != 0 ) answer.times(1. / norm);
        //    }
        //}
        //else {
        //    //todo maybe force reestablishing closest edge? instead of automatically giving up
        //    answer.resize(2);
        //}
    }

    void ElementEdgeContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //return zeros
            answer.resize(2, 6);
            return;
        }

        StructuralElement* elem = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        int edgePos = closestEdge.at(2);

        FloatMatrix N;
        FloatArray cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords;
        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes(edgeNodes, edgePos);

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        elem->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
	    elem->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);


        bool inbetween = computeContactPoint(cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords);

        //all the previous just to compute the contact point...

        FloatArray lcoords;
        elem->giveInterpolation()->global2local(lcoords, cPoint, FEIElementGeometryWrapper(elem));
        elem->computeEdgeNMatrix(N, edgePos, lcoords);

        answer.resize(N.giveNumberOfRows(), N.giveNumberOfColumns() + 2);
        answer.zero();

        FloatMatrix extension(2, 2);
        extension.beUnitMatrix();
        extension.times(-1.);

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

        FloatArray cPoint, projection, edgeNode1Coords, edgeNode2Coords, nodeCoords;
        FloatArray cPointInit, projectionInit, edgeNode1CoordsInit, edgeNode2CoordsInit, nodeCoordsInit;
        IntArray edgeNodes;

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        nodeCoordsInit = node->giveNodeCoordinates();

        StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        element->giveBoundaryEdgeNodes(edgeNodes, closestEdge.at(2));

        element->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
	    element->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);
	
        edgeNode1CoordsInit = element->giveNode(edgeNodes(0))->giveNodeCoordinates();
	    edgeNode2CoordsInit = element->giveNode(edgeNodes(1))->giveNodeCoordinates();

        bool inbetween = computeContactPoint(cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords);
        bool inbetweenInit = computeContactPoint(cPointInit, nodeCoordsInit, edgeNode1CoordsInit, edgeNode2CoordsInit);
        projection.beDifferenceOf(cPoint, nodeCoords);
        projectionInit.beDifferenceOf(cPointInit, nodeCoordsInit);

        //test whether initial and current vector are on different sides of line
        //find whether their cosine is larger than zero
        double cos = projection.dotProduct(projectionInit);
        bool penetrated = cos <= 0;

        double answer = projection.computeNorm();
        if ( penetrated ) answer *= -1;

        return answer;
    }

    void ElementEdgeContactSegment::computeNormalSlope(FloatArray & answer, Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //SIZE ????
            answer.resize(6);
            return;
        }
        int edgePos = closestEdge.at(2);

        StructuralElement* elem = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        NLStructuralElement * nlelem = dynamic_cast<NLStructuralElement*> (elem);
        if ( nlelem == nullptr || nlelem->giveGeometryMode() != 1 ) {
            //there is no geometrical non-linearity
            answer.resize(6);
            return;
        }

        FloatArray cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords;
        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes(edgeNodes, edgePos);

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        elem->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
        elem->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);


        bool inbetween = computeContactPoint(cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords);

        //all the previous just to compute the contact point...

        FloatArray lcoords;
        elem->giveInterpolation()->global2local(lcoords, cPoint, FEIElementGeometryWrapper(elem));

        //the term we are supposed to return is
        // 1/||n|| (n_0 . (F X) + (n/||n||) x (F X) : (n/||n||) x n_0) Bh
        //where
        // n is the deformed normal, n_0 the original normal
        // F is the deformation gradient at contact point
        // Bh is the BH matrix at contact point
        // . represents single contraction, : double contraction, x a dyadic product and X a tensor cross product

        FloatArray n, n_norm, n0; // 2x1
        FloatArray Fv; // 6x1
        double n_size = 0.;
        FloatMatrix Bh;
        FloatMatrix F; //3x3
        FloatMatrix F_cross; //9x9

        nlelem->computeDeformationGradientVector(Fv, lcoords, tStep);

        //retrieve edge normal from element interpolation
        FEInterpolation2d* interpolation = dynamic_cast<FEInterpolation2d*>(elem->giveInterpolation());
        if ( interpolation == nullptr ) {
            OOFEM_ERROR("Non-2D element encountered in ElementEdgeContactSegment");
        }
        interpolation->edgeEvalNormal(n0, edgePos, lcoords, FEIElementGeometryWrapper(elem));
        n = n0;
        transformNormalToDeformedShape(n, nlelem, lcoords, tStep);
        n_size = n.computeNorm();
        n_norm = n;
        n_norm.times(1. / n_size);

        StructuralMaterial::compute_tensor_cross_product_tensor(F_cross, Fv);
        nlelem->computeBHmatrixAt(lcoords, Bh);

        //TODO assemble together. Do tensor sizes agree??
    }

    void ElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s)
    {
        if ( lastEdge.giveSize() == 2 ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(lastEdge.at(1));
            int edgePos = lastEdge.at(2);

            IntArray boundaryNodes;
            element->giveBoundaryEdgeNodes(boundaryNodes, edgePos);
            element->giveBoundaryLocationArray(answer, boundaryNodes, c_s, nullptr);
        } else {
            //if no segment was worked with, returns zeros
            answer.resize(4);
        }

    }

    void ElementEdgeContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s)
    {
        answer.resize(0);
        IntArray edgeloc, boundaryNodes;
        //iterate over all edges and add their locarrays
        for ( int pos = 0; pos < edges.giveSize() / 2; pos++ ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(edges(pos * 2));
            int edgePos = pos * 2 + 1;

            element->giveBoundaryEdgeNodes(boundaryNodes, edges(edgePos));
            element->giveBoundaryLocationArray(edgeloc, boundaryNodes, c_s, nullptr);

            answer.followedBy(edgeloc);
        }
    }

    void ElementEdgeContactSegment::transformNormalToDeformedShape(FloatArray& normal, NLStructuralElement * elem, const FloatArray& lcoords, TimeStep* tStep)
    {
        FloatArray F;
        elem->computeDeformationGradientVector(F, lcoords, tStep);
        //compute cofactor and multiply by normal

        FloatMatrix cofactor;
        StructuralMaterial::compute_2order_tensor_cross_product(cofactor, F, F);
        cofactor.times(0.5);

        normal.beProductOf(cofactor, normal);
    }

    void ElementEdgeContactSegment::giveClosestEdge(IntArray & answer, Node * node, TimeStep * tStep)
    {
        int knownIndex = giveIndexOfKnownNode(node);
        if (updateMode != UM_EveryQuery && knownIndex != -1 && knownClosestEdges.at(knownIndex).giveSize() == 2 ) {
            answer = knownClosestEdges.at(knownIndex);
            lastEdge = answer;
            return;
        }
        //if previous if failed, it means that we don't know the closest edge to this node yet
        answer.resize(2);

        FloatArray cPoint, projection, edgeNode1Coords, edgeNode2Coords, nodeCoords;
        IntArray edgeNodes;
        double answerSize = -1.;

        node->giveUpdatedCoordinates(nodeCoords, tStep);

        //iterate over all edges to find the closest one
        for ( int pos = 0; pos < edges.giveSize() / 2; pos++ ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(edges(pos * 2));
            int edgePos = pos * 2 + 1;

            element->giveBoundaryEdgeNodes(edgeNodes, edges(edgePos));

	        element->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
	        element->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);


            bool inbetween = computeContactPoint(cPoint, nodeCoords, edgeNode1Coords, edgeNode2Coords);
            projection.beDifferenceOf(cPoint, nodeCoords);

            if ( inbetween ) {
                double normalSize = projection.computeNorm();
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
            IntArray answerToStore(answer); //make copy to store to prevent it going out of scope - necessary??
            knownClosestEdges.push_back(answerToStore);
        }
        lastEdge = answer;
    }

    bool ElementEdgeContactSegment::computeContactPoint(FloatArray & answer, const FloatArray & externalPoint, const FloatArray & linePoint1, const FloatArray & linePoint2)
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
        double ae = normalVector.at(1);
        double be = normalVector.at(2);
        double ce = -ae * linePoint1.at(1) - be * linePoint1.at(2);

        //now we have to find the intersecting point
        double divider = ae * bn - an * be;//should only be zero if lines are parallel
        //FloatArray contactPoint(2);
        answer.resize(2);
        answer.at(1) = -(ce*bn - cn * be) / divider;
        answer.at(2) = -(ae*cn - an * ce) / divider;

        //answer.beDifferenceOf(contactPoint, externalPoint);
        //now the contact point is the answer rather than the projection

        double lineLength = abs(linePoint1(0) - linePoint2(0));
        //point lies inbetween line points if it is closer to both than length of line
        return abs(answer(0) - linePoint1(0)) <= lineLength && abs(answer(0) - linePoint2(0)) <= lineLength;
    }

    void ElementEdgeContactSegment::updateYourself(TimeStep *tStep)
        // Updates the receiver at end of step.
    {
        if ( updateMode != UM_Never ) {
            knownNodes.clear();
            knownClosestEdges.clear();
        }
        lastEdge.resize(0);
    }


    void ElementEdgeContactSegment::postInitialize()
    {

        Set* set = this->giveDomain()->giveSet(this->setnum);
        if ( set == nullptr ) OOFEM_ERROR("Contact segment can not find set no. " + setnum);
        this->edges = set->giveBoundaryList();
        if ( edges.giveSize() <= 0 ) OOFEM_WARNING("Contact segment's edge list is empty");

    }



}


