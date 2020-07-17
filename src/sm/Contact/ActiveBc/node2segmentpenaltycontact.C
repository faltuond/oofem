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


#include "sm/Contact/ActiveBc/node2segmentpenaltycontact.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition( Node2SegmentPenaltyContact );


IRResultType Node2SegmentPenaltyContact::initializeFrom( InputRecord *ir )
{
    IRResultType result;
    IR_GIVE_FIELD( ir, this->penalty, _IFT_Node2SegmentPenaltyContact_penalty );
    this->useTangent = ir->hasField( _IFT_Node2SegmentPenaltyContact_useTangent );


    IR_GIVE_FIELD( ir, this->segmentSet, _IFT_Node2SegmentPenaltyContact_segmentSet );
    IR_GIVE_FIELD( ir, this->nodeSet, _IFT_Node2SegmentPenaltyContact_nodeSet );

    IR_GIVE_OPTIONAL_FIELD( ir, this->prescribedNormal, _IFT_Node2SegmentPenaltyContact_prescribedNormal );

    return ActiveBoundaryCondition::initializeFrom( ir );
}


void Node2SegmentPenaltyContact::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale )
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }

    FloatMatrix K;
    IntArray loc, node_loc;

    IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
    /*IntArray dofIdArray = {
          D_u, D_v
        };*/

    //iterate over all pairs of nodes and segments
    for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
        for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {

            Node *node              = this->giveDomain()->giveNode( nodeSet.at( nodePos ) );
            ContactSegment *segment = ( this->giveDomain()->giveContactSegment( segmentSet.at( segmentPos ) ) );

            this->computeTangentFromContact( K, node, segment, tStep );

            //assembling for both node and segment
            node->giveLocationArray( dofIdArray, node_loc, r_s );
            segment->giveLocationArray( dofIdArray, loc, r_s );
            loc.followedBy( node_loc );

            answer.assemble( loc, K );
            //K.printYourself();
        }
    }
}

void Node2SegmentPenaltyContact::assembleVector( FloatArray &answer, TimeStep *tStep,
    CharType type, ValueModeType mode,
    const UnknownNumberingScheme &s, FloatArray *eNorms )
{
    if ( type != ExternalForcesVector ) {
        return;
    }

    IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
    /*IntArray dofIdArray = {
          D_u, D_v
        };*/

    IntArray loc, node_loc;
    FloatArray fext;

    for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
        for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
            Node *node              = this->giveDomain()->giveNode( nodeSet.at( nodePos ) );
            ContactSegment *segment = ( this->giveDomain()->giveContactSegment( segmentSet.at( segmentPos ) ) );

            this->computeExternalForcesFromContact( fext, node, segment, tStep );

            //assembling for both node and segment
            node->giveLocationArray( dofIdArray, node_loc, s );
            segment->giveLocationArray( dofIdArray, loc, s );
            loc.followedBy( node_loc );

            answer.assemble( fext, loc );
        }
    }
}

void Node2SegmentPenaltyContact::giveLocationArrays( std::vector<IntArray> &rows, std::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s )
{
    //returns all possible combinations of dof that can theoretically be triggered by contact
    //of any segment with any node. Room for optimization aplenty...
    IntArray n_loc, s_loc;

    int ncombinations = nodeSet.giveSize() * segmentSet.giveSize();
    rows.resize( ncombinations * 2 );
    cols.resize( ncombinations * 2 );

    IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
    /*IntArray dofIdArray = {
          D_u, D_v
        };*/

    int pos = 0;

    for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); nodePos++ ) {
        for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
            Node *node              = this->giveDomain()->giveNode( nodeSet.at( nodePos ) );
            ContactSegment *segment = ( this->giveDomain()->giveContactSegment( segmentSet.at( segmentPos ) ) );

            node->giveLocationArray( dofIdArray, n_loc, r_s );
            segment->giveLocationArrays( dofIdArray, s_loc, c_s );

            // insert location arrays into the answer arrays
            rows[pos]     = n_loc;
            cols[pos]     = s_loc;
            rows[pos + 1] = s_loc;
            cols[pos + 1] = n_loc;
            pos += 2;
        }
    }
}


void Node2SegmentPenaltyContact::computeGap( double &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    answer = segment->computePenetration( node, tStep );
}


void Node2SegmentPenaltyContact::computeNormalMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{

    FloatArray normal;
    FloatMatrix extendedN, extendedNTranspose;

    if ( prescribedNormal.giveSize() == node->giveNumberOfDofs() ) {
        normal = prescribedNormal;
        /* double gap;
            computeGap(gap, node, segment, tStep);
            if ( gap < 0 ) normal.times(-1.);*/
    } else {
        segment->computeNormal( normal, node, tStep );
    }


    segment->computeExtendedNMatrix( extendedN, node, tStep );
    extendedNTranspose.beTranspositionOf( extendedN );

    //normal should be given just as N^t * n;
    answer.beProductOf( extendedNTranspose, normal );
}

void Node2SegmentPenaltyContact::computeTangentMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    FloatArray tangent;
    FloatMatrix extendedN;

    if ( prescribedNormal.giveSize() == node->giveNumberOfDofs() ) {
        OOFEM_WARNING( "Prescribed normal inapplicable for use with large strains" );
    }
    segment->computeTangent( tangent, node, tStep );

    segment->computeExtendedNMatrix( extendedN, node, tStep );

    //tangent should be given just as N^t * t;
    answer.beTProductOf( extendedN, tangent );
}


void Node2SegmentPenaltyContact::computeTangentFromContact( FloatMatrix &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    double gap;
    this->computeGap( gap, node, segment, tStep );

    //assembling the first part of the tangent
    //considering Nv = N^T * n/||n||
    //and         K1 = p * (Nv*Nv^T)
    // (this part is always assembled as it is the simplest way to know the required dimensions of answer)
    // (if there is no contact, the answer, already with the right dimensions, is later zeroed, see below)
    FloatArray Nv;
    this->computeNormalMatrixAt( Nv, node, segment, tStep );
    answer.beDyadicProductOf( Nv, Nv );
    answer.times( this->penalty );

    if ( gap < 0.0 ) {
        //assembling the second part of the tangent (for large deformations)
        //whether large deformations are to be considered or not is decided
        //by the segment and reflected in the returned normal slope matrix

        //K2 = - p*g*N^T * n_s
        FloatMatrix k2, Ns, Next;
        segment->computeNormalSlope( Ns, node, tStep );
        segment->computeExtendedNMatrix( Next, node, tStep );
        k2.beTProductOf( Next, Ns );
        k2.times( -gap * this->penalty );
        answer.add( k2 );

        /*k2.printYourself();
            answer.printYourself();*/

    } else {
        //zero in the case of no contact occuring
        answer.zero();
    }
}


void Node2SegmentPenaltyContact::computeExternalForcesFromContact( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    double gap;
    this->computeGap( gap, node, segment, tStep );
    this->computeNormalMatrixAt( answer, node, segment, tStep );
    if ( gap < 0.0 ) {
        answer.times( penalty * gap );
    } else {
        answer.times( 0 );
    }
}


} // namespace oofem
