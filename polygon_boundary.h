#ifndef POLYGON_BOUNDARY_H
#define POLYGON_BOUNDARY_H

#include <stdexcept>
#include <vector>
#include <functional>

//! Point Structure for 2D plane
struct PolygonPoint {
    double x, y;
    PolygonPoint(double x = 0, double y = 0) : x(x), y(y) {}
};

//! Type of Boundary
/*! NEUMANN BC shares the algorithm with MIX BC.
 * DIRICHLET is for reducing calculation.
 * I recommend you to declare BCType clearly for acceleration.
 */
enum PolyBoundaryType {
    POLY_DIRICHLET,  // 狄利克雷条件 u = g
    POLY_NEUMANN,    // 诺依曼条件 ∂u/∂n = h
    POLY_MIX
};

//! All Info about Boundary
/*! start end is the position of boundary.
 * type is the type of BC.
 * alpha beta is the parameter in MIX_BC with alpha*u+beta*du/dn=f(x,y).
 * alpha beta can be omitted in other BC.
 * func is the function of BC.
 */
struct PolygonEdge {
    //! start position of boundary
    PolygonPoint start;
    //! end position of boundary
    PolygonPoint end;
    //! type of BC
    PolyBoundaryType type;
    //! factor of u
    /*! can be omitted */
    double alpha;
    //! factor of du/dn
    /*! can be omitted */
    double beta;
    //! function of BC
    std::function<double(double x, double y, double t)> func;

    //! normal declare
    PolygonEdge(){};
    //! declare with no parameter
    PolygonEdge(PolygonPoint start, PolygonPoint end, PolyBoundaryType type,double (*func)(double x, double y)):start(start),
    end(end), type(type), func([func](double x,double y, double t){return (*func)(x,y);}), alpha(0.0), beta(1.0) {
        if (type == POLY_MIX) {
            throw std::invalid_argument("POLY_MIX not implemented");
        }
    };
    //! declare with time
    PolygonEdge(PolygonPoint start, PolygonPoint end, PolyBoundaryType type,double (*func)(double x, double y, double t)):start(start),
    end(end), type(type), func(func), alpha(0.0), beta(1.0) {
        if (type == POLY_MIX) {
            throw std::invalid_argument("POLY_MIX not implemented");
        }
    };
};

//! Polygon Region defined by Polygon Edge
class PolygonRegion {
private:
    //! edges of polygon
    std::vector<PolygonEdge> edges;  // 多边形边界边
    //! Minimum Rectangle Region that covers the Polygon
    /*! Used for acceleration */
    double x_min, x_max, y_min, y_max;  // 包围盒（用于快速判断）

    //! Judge if the point is in the PolygonRegion.
    bool isPointInside(const PolygonPoint& p) const;

    //! Judge if the point is on the PolygonEdge.
    bool isPointOnEdge(const PolygonPoint& p, PolygonEdge& edge_out, double eps) const;

public:

    //! Declare the PolygonRegion with edges.
    PolygonRegion(const std::vector<PolygonEdge>& edges);

    //! Judge if the point is inside.
    bool isInside(double x, double y) const;

    //! Judge if the Point is on the Boundary (NOT contain time).
    bool isOnBoundary(double x, double y, PolyBoundaryType& type, double& value, double& nx, double& ny,
        double& alpha, double& beta, double eps = 1e-4) const;
    //! Judge if the Point is on the Boundary (Contain time).
    bool isOnBoundary(double x, double y, double t, PolyBoundaryType& type, double& value, double& nx, double& ny,
        double& alpha, double& beta, double eps = 1e-4) const;
};

#endif // POLYGON_BOUNDARY_H
