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
		if (this->hasNonLinearGeometry(node, tStep)) {
			//determine normal from tangent as t vector cross product e3
			FloatArray tangent, tangent3D, e3, normal3D;
			computeTangent(tangent, node, tStep);
			tangent3D.resize(3);
			tangent3D.addSubVector(tangent, 1);
			e3.resize(3);
			e3.at(3) = 1.;
			normal3D.beVectorProductOf(tangent3D, e3);
			answer = normal3D;
			answer.resizeWithValues(2);
		}
		else {
			//determine normal from edge
			IntArray closestEdge;
			giveClosestEdge(closestEdge, node, tStep);
			if (closestEdge.giveSize() != 2) {
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
			if (interpolation == nullptr) {
				OOFEM_ERROR("Non-2D element encountered in ElementEdgeContactSegment");
			}

			interpolation->global2local(cPointLocal, cPoint, FEIElementGeometryWrapper(elem));
			interpolation->edgeEvalNormal(normal, edgePos, cPointLocal, FEIElementGeometryWrapper(elem));

			//contact normal is defined as going towards the edge
			//normal.times(-1.);
			answer = normal;
		}
    }

    void ElementEdgeContactSegment::computeTangent( FloatArray &answer, Node *node, TimeStep *tStep )
    {
        IntArray closestEdge;
        giveClosestEdge( closestEdge, node, tStep );
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //return zeros
            answer.resize( 2 );
            return;
        }

		StructuralElement *elem = (StructuralElement *)this->giveDomain()->giveElement( closestEdge.at( 1 ) );
        int edgePos = closestEdge.at( 2 );
        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes( edgeNodes, edgePos );
        FloatArray nodalCoords;
        FloatArray nodalCoords2;

		elem->giveNode( edgeNodes( 0 ) )->giveUpdatedCoordinates( nodalCoords, tStep );
        elem->giveNode( edgeNodes( 1 ) )->giveUpdatedCoordinates( nodalCoords2, tStep );

		nodalCoords.resizeWithValues( 4 ); //why does FloatArray not have followedBy()??
        nodalCoords.at( 3 ) = nodalCoords2.at( 1 );
        nodalCoords.at( 4 ) = nodalCoords2.at( 2 );

        FloatMatrix dNextdXi;
        computeExtendedBMatrix( dNextdXi, node, tStep );
		dNextdXi.resizeWithData(2, 4);

		answer.beProductOf( dNextdXi, nodalCoords );
		answer.times(-1.);
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

    void ElementEdgeContactSegment::computeExtendedBMatrix( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
		//for linear segments, this is always the same
		FloatMatrix answerT;
		answerT = { {1, 0, -1,  0, 0, 0},
				    {0, 1,  0, -1, 0, 0}};
        answer.beTranspositionOf( answerT );
    }

    double ElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
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

    bool ElementEdgeContactSegment::hasNonLinearGeometry( Node *node, TimeStep *tStep )
    {
        IntArray closestEdge;
        giveClosestEdge( closestEdge, node, tStep );
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            return false;
        }

		int edgePos = closestEdge.at( 2 );

        NLStructuralElement *elem = dynamic_cast<NLStructuralElement *> (this->giveDomain()->giveElement( closestEdge.at( 1 ) ));
        return elem != nullptr && elem->giveGeometryMode() == 1;
    }

    void ElementEdgeContactSegment::computeMetricTensor( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
        answer.resize( 1, 1 );
        //the answer is length of segment, which is in fact the size of tangent, squared
		FloatArray tangent;
        computeTangent( tangent, node, tStep );
		double l = tangent.computeNorm();
        answer.at( 1, 1 ) = l * l;
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


