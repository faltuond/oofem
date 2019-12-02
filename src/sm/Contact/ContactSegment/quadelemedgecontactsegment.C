#pragma once

#include "quadelemedgecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(QuadraticElementEdgeContactSegment);

    IRResultType QuadraticElementEdgeContactSegment::initializeFrom(InputRecord* ir) {

        return ContactSegment::initializeFrom(ir);
    }

    void QuadraticElementEdgeContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tstep)
    {
    }

    void QuadraticElementEdgeContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
    }

    double QuadraticElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        return 0.0;
    }

    void QuadraticElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
    }

    void QuadraticElementEdgeContactSegment::updateYourself(TimeStep * tStep)
    {
    }

    void QuadraticElementEdgeContactSegment::postInitialize()
    {
    }

    void QuadraticElementEdgeContactSegment::giveClosestEdge(IntArray & answer, Node * node, TimeStep * tStep)
    {
    }

}