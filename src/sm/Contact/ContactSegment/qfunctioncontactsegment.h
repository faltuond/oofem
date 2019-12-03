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
#include "classfactory.h"
#include "nrfunctioncontactsegment2d.h"

#define _IFT_QFunctionContactSegment_Name "qfunctioncontactsegment"
#define _IFT_QFunctionContactSegment_a "a"
#define _IFT_QFunctionContactSegment_b "b"
#define _IFT_QFunctionContactSegment_c "c"

namespace oofem {
    class QFunctionContactSegment : public NRFunctionContactSegment2D
    {
    public:
        QFunctionContactSegment(int n, Domain *aDomain) : NRFunctionContactSegment2D(n, aDomain) { ; }
        ~QFunctionContactSegment() {};

        IRResultType initializeFrom(InputRecord * ir) override;

        const char *giveClassName() const override { return "Qfunctioncontactsegment"; }
        const char *giveInputRecordName() const override { return _IFT_QFunctionContactSegment_Name; }

    private:
        double a, b, c;

    protected:

        inline double functionValue(const double x) const override { return ((a*x*x) + b*x + c); };
        inline double derivativeValue(const double x) const override { return (2*a*x + b); };
        inline double doubleDerivativeValue(const double x) const override { return 2 * a; };
    };

}