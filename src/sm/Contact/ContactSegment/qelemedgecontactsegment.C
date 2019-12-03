#pragma once

#include "qelemedgecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(QElementEdgeContactSegment);

    IRResultType QElementEdgeContactSegment::initializeFrom(InputRecord* ir) {

        return ContactSegment::initializeFrom(ir);
    }

    void QElementEdgeContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tstep)
    {
    }

    void QElementEdgeContactSegment::computeExtendedNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
    }

    double QElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        return 0.0;
    }

    void QElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
    }

    void QElementEdgeContactSegment::updateYourself(TimeStep * tStep)
    {
    }

    void QElementEdgeContactSegment::postInitialize()
    {
    }

    void QElementEdgeContactSegment::giveClosestEdge(IntArray & answer, Node * node, TimeStep * tStep)
    {
    }

}