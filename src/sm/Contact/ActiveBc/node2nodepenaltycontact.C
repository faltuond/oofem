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


#include "sm/Contact/ActiveBc/node2nodepenaltycontact.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"

#ifdef _OPENMP
 #include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition(Node2NodePenaltyContact);



void
Node2NodePenaltyContact::initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition::initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->penalty, _IFT_Node2NodePenaltyContact_penalty);
    this->useTangent = ir.hasField(_IFT_Node2NodePenaltyContact_useTangent);


    IR_GIVE_FIELD(ir, this->masterSetNumber, _IFT_Node2NodePenaltyContact_masterSet);
    IR_GIVE_FIELD(ir, this->slaveSetNumber, _IFT_Node2NodePenaltyContact_slaveSet);

    IR_GIVE_OPTIONAL_FIELD(ir, this->prescribedNormal, _IFT_Node2NodePenaltyContact_prescribedNormal);
}

void
Node2NodePenaltyContact::postInitialize()
{
    masterSet = domain->giveSet(this->masterSetNumber)->giveNodeList();
    slaveSet  = domain->giveSet(this->slaveSetNumber)->giveNodeList();
}




void
Node2NodePenaltyContact::assemble(SparseMtrx &answer, TimeStep *tStep,
                                  CharType type, const UnknownNumberingScheme &r_s,
                                  const UnknownNumberingScheme &c_s, double scale,
                                  void *lock)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }


    FloatMatrix K;
    IntArray loc, c_loc;

    IntArray dofIdArray = domain->giveDefaultNodeDofIDArry();

    if ( masterSet.giveSize() != slaveSet.giveSize() ) {
        OOFEM_ERROR("Number of master nodes does not match number of slave nodes")
    }

    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode( masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode( slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, loc, r_s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);
        loc.followedBy(c_loc);

        this->computeTangentFromContact(K, masterNode, slaveNode, tStep);
#ifdef _OPENMP
        if ( lock ) {
            omp_set_lock( static_cast< omp_lock_t * >( lock ) );
        }
#endif
        answer.assemble(loc, K);
#ifdef _OPENMP
        if ( lock ) {
            omp_unset_lock( static_cast< omp_lock_t * >( lock ) );
        }
#endif
    }
}

void
Node2NodePenaltyContact::assembleVector(FloatArray &answer, TimeStep *tStep,
                                        CharType type, ValueModeType mode,
                                        const UnknownNumberingScheme &s, FloatArray *eNorms,
                                        void *lock)
{
    if ( type != ExternalForcesVector ) {
        return;
    }

    //IntArray dofIdArray = {D_u, D_v, D_w};
    IntArray dofIdArray = domain->giveDefaultNodeDofIDArry();


    IntArray loc, c_loc;
    FloatArray fext;

    if ( masterSet.giveSize() != slaveSet.giveSize() ) {
        OOFEM_ERROR("Number of master nodes does not match number of slave nodes");
    }

    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode( masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode( slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, loc, s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, s);
        this->computeExternalForcesFromContact(fext, masterNode, slaveNode, tStep);
        loc.followedBy(c_loc);
#ifdef _OPENMP
        if ( lock ) {
            omp_set_lock( static_cast< omp_lock_t * >( lock ) );
        }
#endif
        answer.assemble(fext, loc);
#ifdef _OPENMP
        if ( lock ) {
            omp_unset_lock( static_cast< omp_lock_t * >( lock ) );
        }
#endif
    }
}



void Node2NodePenaltyContact::giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray r_loc, c_loc;
    rows.resize(masterSet.giveSize() * 2);
    cols.resize(masterSet.giveSize() * 2);
    /*IntArray dofIdArray = {
        D_u, D_v
    };*/
    IntArray dofIdArray = domain->giveDefaultNodeDofIDArry();

    int arraypos = 0;
    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode( masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode( slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, r_loc, r_s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);

        // column block
        rows [ arraypos ] = r_loc;
        cols [ arraypos ] = c_loc;
        rows [ arraypos + 1 ] = c_loc;
        cols [ arraypos + 1 ] = r_loc;

        //arraypos++;
        arraypos += 2;
    }
}


void Node2NodePenaltyContact::computeTangentFromContact(FloatMatrix &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    double gap;
    FloatArray Nv;
    this->computeGap(gap, masterNode, slaveNode, tStep);
    this->computeNvMatrixAt(Nv, masterNode, slaveNode, tStep);
    answer.beDyadicProductOf(Nv, Nv);
    answer.times(this->penalty);
    if ( gap > 0.0 ) {
        answer.zero();
    }
}

void Node2NodePenaltyContact::computeGap(double &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    IntArray dofIdArray = domain->giveDefaultNodeDofIDArry();

    FloatArray uS, uM;
    auto xs = slaveNode->giveCoordinates();
    auto xm = masterNode->giveCoordinates();
    FloatArray normal;
    if ( prescribedNormal.giveSize() == masterNode->giveNumberOfDofs() ) {
        normal = prescribedNormal;
    } else   {
        normal = xs - xm;
    }
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR( "Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.",
                     masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else   {
        normal.times(1.0 / norm);
    }

    slaveNode->giveUnknownVector(uS, dofIdArray, VM_Total, tStep, true);
    masterNode->giveUnknownVector(uM, dofIdArray, VM_Total, tStep, true);
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs - xm;
    answer = dx.dotProduct(normal);
}


void Node2NodePenaltyContact::computeNvMatrixAt(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *TimeStep)
{
    //int dimension = masterNode->giveNumberOfDofs();

    FloatArray normal;
    if ( prescribedNormal.giveSize() == masterNode->giveNumberOfDofs() ) {
        normal = prescribedNormal;
    } else   {
        const auto &xs = slaveNode->giveCoordinates();
        const auto &xm = masterNode->giveCoordinates();
        normal = xs - xm;
    }
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR( "Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.",
                     masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else   {
        normal.times(1.0 / norm);
    }
    // The normal is not updated for node2node which is for small deformations only
    // C = {n -n}
    //answer = {
    //    normal.at(1), normal.at(2),
    //    -normal.at(1), -normal.at(2)
    //};

    FloatArray negativeNormal;
    negativeNormal = normal;
    negativeNormal.times(-1.);
    answer.resize(2 * normal.giveSize());
    answer.addSubVector(normal, 1);
    answer.addSubVector(negativeNormal, normal.giveSize() + 1);
}


void Node2NodePenaltyContact::computeExternalForcesFromContact(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    double gap;
    FloatArray fext;
    this->computeGap(gap, masterNode, slaveNode, tStep);
    this->computeNvMatrixAt(fext, masterNode, slaveNode, tStep);
    if ( gap < 0.0 ) {
        fext.times(penalty * gap);
    } else   {
        fext.times(0);
    }
    answer = fext;
}
} // namespace oofem
