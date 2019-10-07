/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "node2segmentpenaltycontact.h"

namespace oofem {

    Node2SegmentPenaltyContact::Node2SegmentPenaltyContact(int n, Domain * d) : ActiveBoundaryCondition(n, d)
    {
    }
    
    Node2SegmentPenaltyContact::~Node2SegmentPenaltyContact()
    {
    }
    
    IRResultType Node2SegmentPenaltyContact::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        IR_GIVE_FIELD(ir, this->penalty, _IFT_Node2SegmentPenaltyContact_penalty);
        this->useTangent = ir->hasField(_IFT_Node2SegmentPenaltyContact_useTangent);


        IR_GIVE_FIELD(ir, this->segmentSet, _IFT_Node2SegmentPenaltyContact_segmentSet);
        IR_GIVE_FIELD(ir, this->nodeSet, _IFT_Node2SegmentPenaltyContact_nodeSet);
        
        return ActiveBoundaryCondition::initializeFrom(ir);
    }

    void Node2SegmentPenaltyContact::assemble(SparseMtrx & answer, TimeStep * tStep, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s, double scale)
    {
        if ( !this->useTangent || type != TangentStiffnessMatrix ) {
            return;
        }
        
        FloatMatrix K;
        IntArray loc;

        IntArray dofIdArray = {
            D_u, D_v
        };

        //iterate over all pairs of nodes and segments
        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {

                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));

                //check whether element is a valid contact segment
                StructuralElement* element = dynamic_cast<StructuralElement*>(this->giveDomain()->giveElement(segmentSet.at(segmentPos)));
                Node2SegmentInterface* segment = dynamic_cast<Node2SegmentInterface*>(element);

                if ( segment == nullptr ) {
                    OOFEM_ERROR("A specified contact element is not an instance of Node2SegmentInterface");
                    return;
                }

                //getting locarrays just from node, not from element. Problem?
                //probable course of action - element shall assemble its part by itself
                //question - how does element know the penalty value??
                node->giveLocationArray(dofIdArray, loc, r_s);

                this->computeTangentFromContact(K, node, segment, tStep);
                answer.assemble(loc, K);
            }
        }
    }
    void Node2SegmentPenaltyContact::assembleVector(FloatArray & answer, TimeStep * tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme & s, FloatArray * eNorms)
    {
        if ( type != ExternalForcesVector ) {
            return;
        }

        //IntArray dofIdArray = {D_u, D_v, D_w};
        IntArray dofIdArray = {
            D_u, D_v
        };

        IntArray loc;
        FloatArray fext;

        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));

                //check whether element is a valid contact segment
                StructuralElement* element = dynamic_cast<StructuralElement*>(this->giveDomain()->giveElement(segmentSet.at(segmentPos)));
                Node2SegmentInterface* segment = dynamic_cast<Node2SegmentInterface*>(element);

                if ( segment == nullptr ) {
                    OOFEM_ERROR("A specified contact element is not an instance of Node2SegmentInterface");
                    return;
                }

                //again - loc just from node??
                node->giveLocationArray(dofIdArray, loc, s);
                this->computeExternalForcesFromContact(fext, node, segment, tStep);
                answer.assemble(fext, loc);
            }
                       
        }
    }
    void Node2SegmentPenaltyContact::computeTangentFromContact(FloatMatrix & answer, Node * node, Node2SegmentInterface * segment, TimeStep * tStep)
    {
        double gap;
        FloatArray Nv;
        this->computeGap(gap, node, segment, tStep);
        this->computeNormalMatrixAt(Nv, node, segment, tStep);
        answer.beDyadicProductOf(Nv, Nv);
        answer.times(this->penalty);
        if ( gap > 0.0 ) answer.zero();
    }
    void Node2SegmentPenaltyContact::computeGap(double & answer, Node *node, Node2SegmentInterface *segment, TimeStep * tStep)
    {
        answer = segment->computePenetration(node);
    }
    void Node2SegmentPenaltyContact::computeNormalMatrixAt(FloatArray & answer, Node * node, Node2SegmentInterface * segment, TimeStep * TimeStep)
    {
        //this implementation assumes that the normal computed by Node2SegmentInterface::computeProjection()
        //is the normal in the undeformed state, and that this is the normal to be used all the time
        //i. e. small deformations are assumed

        //if Node2SegmentInterface::computeProjection() would compute the normal in deformed state, problems could
        //arise when the normal is normalized (divide by (a number close to) zero)

        FloatArray normal;
        segment->computeProjection(normal, node);
        double norm = normal.computeNorm();
        normal.times(1.0 / norm);
        
        answer = {
            normal.at(1), normal.at(2),
            -normal.at(1), -normal.at(2)
        };
    }
    void Node2SegmentPenaltyContact::computeExternalForcesFromContact(FloatArray & answer, Node * node, Node2SegmentInterface * segment, TimeStep * tStep)
    {
        double gap;
        this->computeGap(gap, node, segment, tStep);
        this->computeNormalMatrixAt(answer, node, segment, tStep);
        if ( gap < 0.0 ) {
            answer.times(penalty * gap);
        }
        else {
            answer.times(0);
        }
    }
    void Node2SegmentPenaltyContact::giveLocationArrays(std::vector<IntArray>& rows, std::vector<IntArray>& cols, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s)
    {
        //returns all possible combinations of dof that can theoretically be triggered by contact
        //of any segment with any node. Room for optimization aplenty...
        IntArray n_loc, s_loc;

        int ncombinations = nodeSet.giveSize() * segmentSet.giveSize();
        rows.resize(ncombinations);
        cols.resize(ncombinations);
        IntArray dofIdArray = {
            D_u, D_v
        };

        int pos = 0;

        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); nodePos++ ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));
                StructuralElement* element = dynamic_cast<StructuralElement*>(this->giveDomain()->giveElement(segmentSet.at(segmentPos)));
                Node2SegmentInterface* segment = dynamic_cast<Node2SegmentInterface*>(element);
                if ( segment == nullptr ) {
                    OOFEM_ERROR("A specified contact element is not an instance of Node2SegmentInterface");
                    return;
                }

                node->giveLocationArray(dofIdArray, n_loc, r_s);
                segment->giveLocationArray(dofIdArray, s_loc, c_s);

                // insert location arrays into the answer fields
                rows[pos] = n_loc;
                cols[pos] = s_loc;
                pos++;
            }
        }
    }
}