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

#pragma once

#include "activebc.h"
#include "node.h"
#include "floatmatrix.h"
#include "sm/Contact/ContactSegment/ContactSegment.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"
#include "sm/Elements/structuralelement.h"

#define _IFT_Node2SegmentPenaltyContact_Name "n2spenaltycontact"
#define _IFT_Node2SegmentPenaltyContact_penalty "penalty"
#define _IFT_Node2SegmentPenaltyContact_useTangent "usetangent"

#define _IFT_Node2SegmentPenaltyContact_segmentSet "segmentset"
#define _IFT_Node2SegmentPenaltyContact_nodeSet "nodeset"

namespace oofem {
    class Node2SegmentPenaltyContact : ActiveBoundaryCondition
    {
    public:
        Node2SegmentPenaltyContact(int n, Domain *d);
        ~Node2SegmentPenaltyContact();

        bool useTangent; ///< Determines if tangent should be used.
        double penalty;
        IntArray nodeSet;
        IntArray segmentSet;

        virtual IRResultType initializeFrom(InputRecord *ir);

        virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
            CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

        virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
            CharType type, ValueModeType mode,
            const UnknownNumberingScheme &s, FloatArray *eNorms = NULL) override;


        virtual const char *giveClassName() const { return "Node2SegmentPenaltyContact"; }
        virtual const char *giveInputRecordName() const { return _IFT_Node2SegmentPenaltyContact_Name; }


        void computeTangentFromContact(FloatMatrix &answer, Node *node, ContactSegment *segment, TimeStep *tStep);
        void computeGap(double &answer, Node *node, ContactSegment *segment, TimeStep *tStep);

        void computeNormalMatrixAt(FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *TimeStep);


        void computeExternalForcesFromContact(FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep);

        void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    };

}