#include "sparse_matrix.h"
#include "PoissonEQ.h"
#include <numeric>
#include <chrono>
#include <cmath>
#include "ParabolaPDE.h"
#include "HyperbolicEQ.h"

double delta(double x, double y, double x0, double y0, double tol) {
    double delX = x-x0;
    double delY = y-y0;
    if (delX*delX + delY*delY < tol) {
        return 1;
    }else {
        return 0;
    }
}

double initCondition(double x, double y) {
    return 0.0;
}

double sinrr(double beta, double x) {
    if (abs(x)<1e-5) {
        return beta;
    }
    return sin(beta*x)/x;
}

double sourceTerm(double x, double y, double t) {
    double r = sqrt(x*x + y*y);
    return exp(-5.0*r*r)*(50*cos(2*r)*cos(50*t)-100*cos(50*r));
}

double testFunc(double x, double y, double t) {
    return 1/(4*std::numbers::pi*(t+0.1))*exp(-((x-1)*(x-1)+(y-1)*(y-1))/(4*(t+0.1)));
}

double testInitCondition(double x, double y) {
    return testFunc(x, y, 0.0);
}

double testSourceTerm(double x, double y, double t) {
    return 1/(4*std::numbers::pi*(t+0.1)*(t+0.1))*delta(x,y,1.0,1.0,0.005);
}

int main() {
    solve_hyperbolic_lax_wendroff();
    return 0;
}