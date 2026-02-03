#include "polygon_boundary.h"
#include <cmath>
#include <algorithm>

// 射线法判断点是否在多边形内
bool PolygonRegion::isPointInside(const PolygonPoint& p) const {
    int count = 0;
    for (const auto& edge : edges) {
        const auto& s = edge.start;
        const auto& e = edge.end;

        // 射线与边相交检测
        if (((s.y > p.y) != (e.y > p.y)) &&
            (p.x < s.x + (e.x - s.x) * (p.y - s.y) / (e.y - s.y))) {
            count++;
        }
    }
    return (count % 2) == 1; // 奇数则在内部
}

// 计算点到线段的距离，判断是否在边界上
bool PolygonRegion::isPointOnEdge(const PolygonPoint& p, PolygonEdge& edge_out, double eps) const {
    for (const auto& edge : edges) {
        const auto& s = edge.start;
        const auto& e = edge.end;

        // 计算投影参数
        double dx = e.x - s.x;
        double dy = e.y - s.y;
        double t = ((p.x - s.x) * dx + (p.y - s.y) * dy) / (dx*dx + dy*dy);
        t = std::max(0.0, std::min(1.0, t)); // 限制在[0,1]

        // 计算投影点
        double proj_x = s.x + t * dx;
        double proj_y = s.y + t * dy;

        // 距离小于容差则认为在边界上
        double dist = std::hypot(p.x - proj_x, p.y - proj_y);
        if (dist < eps) {
            edge_out = edge;
            return true;
        }
    }
    return false;
}

PolygonRegion::PolygonRegion(const std::vector<PolygonEdge>& edges) : edges(edges) {
    // 计算包围盒（加速判断）
    x_min = edges[0].start.x;
    x_max = edges[0].start.x;
    y_min = edges[0].start.y;
    y_max = edges[0].start.y;
    for (const auto& edge : edges) {
        x_min = std::min(x_min, std::min(edge.start.x, edge.end.x));
        x_max = std::max(x_max, std::max(edge.start.x, edge.end.x));
        y_min = std::min(y_min, std::min(edge.start.y, edge.end.y));
        y_max = std::max(y_max, std::max(edge.start.y, edge.end.y));
    }
}

bool PolygonRegion::isInside(double x, double y) const {
    // 先通过包围盒快速排除
    if (x < x_min || x > x_max || y < y_min || y > y_max) return false;
    return isPointInside(PolygonPoint(x, y));
}

bool PolygonRegion::isOnBoundary(double x, double y, PolyBoundaryType& type, double& value, double& nx, double& ny,
    double& alpha, double& beta, double eps) const {
    PolygonEdge edge;
    if (isPointOnEdge(PolygonPoint(x, y), edge, eps)) {
        type = edge.type;
        value = edge.func(x, y, 0);
        alpha = edge.alpha;
        beta = edge.beta;

        double lx = edge.end.x - edge.start.x;
        double ly = edge.end.y - edge.start.y;
        nx = ly / sqrt(lx*lx + ly*ly);
        ny = -lx / sqrt(lx*lx + ly*ly);
        return true;
    }
    return false;
}

bool PolygonRegion::isOnBoundary(double x, double y, double t, PolyBoundaryType& type, double& value, double& nx, double& ny,
    double& alpha, double& beta, double eps) const {
    PolygonEdge edge;
    if (isPointOnEdge(PolygonPoint(x, y), edge, eps)) {
        type = edge.type;
        value = edge.func(x, y, t);
        alpha = edge.alpha;
        beta = edge.beta;

        double lx = edge.end.x - edge.start.x;
        double ly = edge.end.y - edge.start.y;
        nx = ly / sqrt(lx*lx + ly*ly);
        ny = -lx / sqrt(lx*lx + ly*ly);
        return true;
    }
    return false;
}
