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
#include "contactsegment.h"
#include "node.h"

#define _IFT_ElementEdgeContactSegment_Name "ElementEdgeContactSegment"
#define _IFT_ElementEdgeContactSegment_edgeSet "edgeset"

namespace oofem {
    class ElementEdgeContactSegment : public ContactSegment
    {
    public:
        ElementEdgeContactSegment(int n, Domain *aDomain) : ContactSegment(n, aDomain) { ; }
        ~ElementEdgeContactSegment() {};

        IRResultType initializeFrom(InputRecord * ir);

        //returns normalized n, which is an normal vector of contact
        virtual void computeNormal(FloatArray& answer, const Node * node) override;

        //returns an extended N (aka A) matrix, integrated at point of contact of given node
        virtual void computeExtendedNMatrix(FloatMatrix& answer, const Node* node) override;

        //computes the penetration of node given 
        virtual double computePenetration(const Node * node) override;

        virtual void giveLocationArray(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) override;
    };
}