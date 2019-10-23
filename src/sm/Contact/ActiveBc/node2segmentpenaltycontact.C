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
        IntArray loc, seg_loc;

        IntArray dofIdArray = {
            D_u, D_v
        };

        //iterate over all pairs of nodes and segments
        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {

                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));
                ContactSegment* segment = (this->giveDomain()->giveContactSegment(segmentSet.at(segmentPos)));

                //assembling for both node and segment
                node->giveLocationArray(dofIdArray, loc, r_s);
                segment->giveLocationArray(dofIdArray, seg_loc, r_s);
                loc.followedBy(seg_loc);

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

        IntArray loc, seg_loc;
        FloatArray fext;

        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));
                ContactSegment* segment = (this->giveDomain()->giveContactSegment(segmentSet.at(segmentPos)));

                //assemble into node and segment loc
                node->giveLocationArray(dofIdArray, loc, s);
                segment->giveLocationArray(dofIdArray, seg_loc, s);
                loc.followedBy(seg_loc);

                this->computeExternalForcesFromContact(fext, node, segment, tStep);
                answer.assemble(fext, loc);
            }
                       
        }
    }
    void Node2SegmentPenaltyContact::computeTangentFromContact(FloatMatrix & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
        //considering that the tangent is given as
        //int across seg (N^T * (n x n) * N), which is equivalent to
        //int across seg ((N^T * n) * (N^T * n)^T)
        //it is sufficient just to dyadically multiply the "normals"

        double gap;
        FloatArray Nv;
        this->computeGap(gap, node, segment, tStep);
        this->computeNormalMatrixAt(Nv, node, segment, tStep);
        answer.beDyadicProductOf(Nv, Nv);
        answer.times(this->penalty);
        if ( gap > 0.0 ) answer.zero();
    }
    void Node2SegmentPenaltyContact::computeGap(double & answer, Node *node, ContactSegment *segment, TimeStep * tStep)
    {
        answer = segment->computePenetration(node);
    }
    void Node2SegmentPenaltyContact::computeNormalMatrixAt(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * TimeStep)
    {
        //computeNormal is expected to return an integrated term
        // int across seg (N^T), where N = element Nmatrix (extended by zeros for the node)

        FloatArray normal;
        FloatMatrix extendedN, extendedNTranspose;
        segment->computeNormal(normal, node);
        segment->computeExtendedNMatrix(extendedN, node);
        extendedNTranspose.beTranspositionOf(extendedN);
        
        //normal should be given just as N^t * n;
        answer.beProductOf(extendedN, normal);
    }
    void Node2SegmentPenaltyContact::computeExternalForcesFromContact(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
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
                ContactSegment* segment = (this->giveDomain()->giveContactSegment(segmentSet.at(segmentPos)));

                node->giveLocationArray(dofIdArray, n_loc, r_s);
                segment->giveLocationArray(dofIdArray, s_loc, c_s);

                // insert location arrays into the answer arrays
                rows[pos] = n_loc;
                cols[pos] = s_loc;
                pos++;
            }
        }
    }
}