#ifndef PARABOLA_CN_H
#define PARABOLA_CN_H

#include "PoissonEQ.h"
#include <vector>
#include <functional>
#include <string>
#include <iostream>

//!  A Parabola Equation in Rectangle Area in CN scheme.
/*!
 * This is a PoissonEQ Class. With these data saved:
 *   X-direction Node Nx
 *   Y-direction Node Ny
 *   X-direction Interval hx
 *   Y-direction Interval hy
 *   Time Interval dt
 *   Conduction factor alpha
 *   A PolyPoissonEQ Solver
 *   PDEMatrix
*/
class ParabolaCN {
private:
    PolyPoissonEQ ellipse_solver;
    SparseMatrix PDEMatrix;

    // CN Scheme parameter
    double alpha;         // Diffusion factor
    double dt;            // Time interval
    int total_steps;      // t steps
    int Nx, Ny;           // Net number
    double hx, hy;        // Space interval



public:
    //! Calculate RHS of CN Form. Need result of the last time step, source function and time range
    std::vector<double> makeRHS(const std::vector<double>& u_prev,
                                  std::function<double(double, double, double)> sourceFunc,
                                  double t_prev, double t_curr);
    //! Save Result
    void saveStepResult(const std::vector<double>& u, int step, std::string filename) const;
    //! Construct Parabola Equation
    ParabolaCN(const PolyPoissonEQ& ellipse_solver,
               double alpha, double dt, int total_steps);
    //! Make LHS of Parabola PDE(NOT contain time)
    void makePDE();
    //! Set Initial Condition. Outer part will be forced to 0.
    std::vector<double> setInitialCondition(std::function<double(double, double)> initFunc);
    //! Iteration solving. Result in all time saved in filename.
    void solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc, std::string filename);
    //! Iteration solving. Result not saved. The last result is returned.
    std::vector<double> solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc);
};


//!  A Parabola Equation in Rectangle Area in explicit scheme.
/*!
 * This is a PoissonEQ Class. With these data saved:
 *   X-direction Node Nx
 *   Y-direction Node Ny
 *   X-direction Interval hx
 *   Y-direction Interval hy
 *   Time Interval dt
 *   Conduction factor alpha
 *   A PolyPoissonEQ Solver
 *   PDEMatrix
*/
class ParabolaExplicit {
private:
    PolyPoissonEQ ellipse_solver;
    SparseMatrix PDEMatrix;

    // Explicit Scheme parameter
    double alpha;         // Diffusion factor
    double dt;            // Time interval
    int total_steps;      // t steps
    int Nx, Ny;           // Net number
    double hx, hy;        // Space interval



public:
    //! Calculate RHS of CN Form. Need result of the last time step, source function and time range
    std::vector<double> makeRHS(const std::vector<double>& u_prev,
                                  std::function<double(double, double, double)> sourceFunc,
                                  double t_curr);
    //! Save Result
    void saveStepResult(const std::vector<double>& u, int step, std::string filename) const;
    //! Construct Parabola Equation
    ParabolaExplicit(const PolyPoissonEQ& ellipse_solver,
               double alpha, double dt, int total_steps);
    //! Make LHS of Parabola PDE(NOT contain time)
    void makePDE();
    //! Set Initial Condition. Outer part will be forced to 0.
    std::vector<double> setInitialCondition(std::function<double(double, double)> initFunc);
    //! Iteration solving. Result in all time saved in filename.
    void solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc, std::string filename);
    //! Iteration solving. Result not saved. The last result is returned.
    std::vector<double> solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc);
};

#endif // PARABOLA_CN_H