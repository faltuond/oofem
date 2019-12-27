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
#include "femcmpnn.h"
#include "node.h"
#include "inputrecord.h"

#define _IFT_ContactSegment_normMode "normmode"

namespace oofem {
    class ContactSegment : public FEMComponent
    {
    public:
      ContactSegment(int n, Domain *aDomain) : FEMComponent(n, aDomain){;}
      ~ContactSegment() {};

      virtual IRResultType initializeFrom(InputRecord* ir) override {
          IRResultType result;

          normmode = NM_Never;
          int normmodeint = 0;
          IR_GIVE_OPTIONAL_FIELD(ir, normmodeint, _IFT_ContactSegment_normMode);
          if ( result == IRRT_OK ) {
              normmode = (NormalizationMode)normmodeint;
              if ( normmodeint < 0 || normmodeint > 2 ) OOFEM_ERROR("Contact segment normalization mode can be only 0, 1 or 2");
          }

          return FEMComponent::initializeFrom(ir);
      };

      //returns normalized n, which is an normal vector of contact
      virtual void computeNormal(FloatArray& answer, Node * node, TimeStep* tstep) = 0;

      //returns an extended N (aka A) matrix, integrated at point of contact of given node
      virtual void computeExtendedNMatrix(FloatMatrix& answer, Node* node, TimeStep * tStep) = 0;

      //computes the penetration of node given 
      virtual double computePenetration(Node * node, TimeStep * tStep) = 0;

      virtual void giveLocationArray(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) = 0;

      virtual void giveLocationArrays(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) = 0;

      virtual void updateYourself(TimeStep * tStep) { ; };

      virtual void postInitialize(){;}

    protected:

        typedef enum NormalizationMode { NM_Never = 0, NM_IfNotSmall = 1, NM_Always = 2 };
        NormalizationMode normmode;
    };
}

