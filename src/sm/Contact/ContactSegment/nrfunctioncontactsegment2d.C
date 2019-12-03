#include "nrfunctioncontactsegment2d.h"

namespace oofem {

    void NRFunctionContactSegment2D::computeDistanceVector(FloatArray & answer, const FloatArray & nodeCoords)
    {
        if ( nodeCoords.giveSize() != 2 ) {
            OOFEM_ERROR("An incompatible coordinate size (%i) encountered. Algorithm is only for 2D functions.");
        }

        double x_node = nodeCoords.at(1);
        double y_node = nodeCoords.at(2);

        //minimizing the distance function d(x) = sqrt((x - x_node)^2 + (f(x) - y_node)^2)
        //omit the unnecessary root and the x_node^2 and y_node^2 to get
        //d(x) = x^2 - 2*x*x_node + f(x)^2 - 2*f(x)*y_node
        //dd(x)/dx = 2*x - 2*x_node + 2* f(x) * df(x)/dx - 2*y_node*df(x)/dx
        //d^2d(x)/dx^2 = 2 + 2*f(x)*(d^2f(x)/dx^2) + 2*(df(x)/dx)^2 - 2*y_node*(d^2f(x)/dx^2)

        //Newton-Rhapson, inital guess x_node
        double x_c = x_node;
        double g, h;
        double k = 0; //iterator
        double maxiter = NRFunctionContact_Maxiter, tol = NRFunctionContact_Tolerance;
        while ( k < maxiter ) {
            g = 2. * x_c - 2. * x_node + 2. * functionValue(x_c)*derivativeValue(x_c) - 2. * y_node*derivativeValue(x_c);
            if ( g <= tol ) break;
            h = 2. + 2.*functionValue(x_c)*doubleDerivativeValue(x_c) + 2.*derivativeValue(x_c)*derivativeValue(x_c) - 2.*y_node*doubleDerivativeValue(x_c);
            x_c += -g / h;
            k++;
        }

        if ( k >= maxiter ) OOFEM_WARNING("Searching for contact with analytical function: Newton-Rhapson method did not converge in %i iterations. Continuing.", k);

        //we found the contact point
        FloatArray contactPoint(x_c, functionValue(x_c));
        answer.beDifferenceOf(contactPoint, nodeCoords);
    }
}