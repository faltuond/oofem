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
    }
    void Node2SegmentPenaltyContact::assembleVector(FloatArray & answer, TimeStep * tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme & s, FloatArray * eNorms)
    {
    }
    void Node2SegmentPenaltyContact::computeTangentFromContact(FloatMatrix & answer, Node * node, Node2SegmentInterface * segment, TimeStep * tStep)
    {
    }
    void Node2SegmentPenaltyContact::computeGap(double & answer, Node * masterNode, Node * slaveNode, TimeStep * tStep)
    {
    }
    void Node2SegmentPenaltyContact::computeNormalMatrixAt(FloatArray & answer, Node * node, Node2SegmentInterface * segment, TimeStep * TimeStep)
    {
    }
    void Node2SegmentPenaltyContact::computeExternalForcesFromContact(FloatArray & answer, Node * node, Node2SegmentInterface * segment, TimeStep * tStep)
    {
    }
    void Node2SegmentPenaltyContact::giveLocationArrays(std::vector<IntArray>& rows, std::vector<IntArray>& cols, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s)
    {
    }
}