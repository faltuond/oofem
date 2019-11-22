#include "node2segmentlagrangemultipliercontact.h"

namespace oofem {
    REGISTER_BoundaryCondition(Node2SegmentLagrangianMultiplierContact);



    IRResultType Node2SegmentLagrangianMultiplierContact::initializeFrom(InputRecord * ir)
    {
        return IRResultType();
    }

    void Node2SegmentLagrangianMultiplierContact::assemble(SparseMtrx & answer, TimeStep * tStep, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s, double scale)
    {
    }

    void Node2SegmentLagrangianMultiplierContact::assembleVector(FloatArray & answer, TimeStep * tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme & s, FloatArray * eNorms)
    {
    }

    void Node2SegmentLagrangianMultiplierContact::giveLocationArrays(std::vector<IntArray>& rows, std::vector<IntArray>& cols, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s)
    {
    }

    void Node2SegmentLagrangianMultiplierContact::giveLagrangianMultiplierLocationArray(const UnknownNumberingScheme & r_s, std::vector<IntArray>& answer)
    {
    }

    double Node2SegmentLagrangianMultiplierContact::computeTangentFromContact(FloatMatrix & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
        return 0.0;
    }

    void Node2SegmentLagrangianMultiplierContact::computeGap(double & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
    }

    void Node2SegmentLagrangianMultiplierContact::computeNormalMatrixAt(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * TimeStep)
    {
    }

    void Node2SegmentLagrangianMultiplierContact::computeExternalForcesFromContact(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
    }

}//end namespace oofem