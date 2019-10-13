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

#include "node2segmentinterface.h"
#include "Elements/Bars/truss2d.h"

#define _IFT_Truss2DNode2SegmentContactElement_Name "truss2dcontactsegment"

namespace oofem {

    class Truss2DNode2SegmentContactElement : Truss2d, Node2SegmentInterface
    {
    public:
        Truss2DNode2SegmentContactElement(int n, Domain* d);
        ~Truss2DNode2SegmentContactElement();

        void computeNormalTerm(FloatArray& answer, const Node* node) override;
        double computePenetration(const Node* node) override;

        const char *giveInputRecordName() const override { return _IFT_Truss2DNode2SegmentContactElement_Name; }
        const char *giveClassName() const override { return "Truss2DNode2SegmentContactElement"; }
    };

}