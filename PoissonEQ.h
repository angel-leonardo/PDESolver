//
// Created by 23923 on 2025/10/5.
//

#ifndef WEEK1_POISSONEQ_H
#define WEEK1_POISSONEQ_H
#include <vector>
#include <fstream>
#include "sparse_matrix.h"
#include "polygon_boundary.h"


//!  A Poisson Equation in Rectangle Area.
/*!
 * This is a PoissonEQ Class. With these data saved:
 *   X-direction Node Nx
 *   Y-direction Node Ny
 *   X-direction Interval hx
 *   Y-direction Interval hy
 *   Differential Form of the EQ
 *   Nonlinear part of the EQ containing Boundary Condition
*/
class PolyPoissonEQ {
private:

    //! Rows of Matrix
    int Nx;
    //! Columns of Matrix
    int Ny;
    //! X Interval
    double hx;
    //! Y Interval
    double hy;
    //! Toleration of Boundary Judgement
    double Boundary_tol;

    //! List of Point to save if the node is active.
    std::vector<bool> isActiveNode;
    //! PolyRegion of EQ
    PolygonRegion* polyRegion;


    //! Nonlinear part of the EQ containing Boundary Condition
    std::vector<double> PDEConstant;
    //! Differential Form of the EQ
    SparseMatrix PDEMatrix;


    public:
    //! Define the RecPoisson
    /*! Need the Information of Grid */
    PolyPoissonEQ(int Nx, int Ny, double hx, double hy);
    //! Construct Poisson Equation in PolygonRegion.
    PolyPoissonEQ(int Nx, int Ny, double hx, double hy, const PolygonRegion& region, double Boundary_tol=1e-4);

    //! Get data in private space
    int getNx() const{return Nx;}
    int getNy() const{return Ny;}
    double gethx() const{return hx;}
    double gethy() const{return hy;}
    double getBoundaryTol() const{return Boundary_tol;}
    std::vector<bool> getActiveNode() const{return isActiveNode;}
    PolygonRegion* getPolyRegion() const{return polyRegion;}

    //! Translate Index to Physical Position X
    double Index2PositionX(int index);
    //! Translate Index to Physical Position Y
    double Index2PositionY(int index);
    //! Translate Index to Node Position X
    int Index2NodeX(int index);
    //! Translate Index to Node Position Y
    int Index2NodeY(int index);


    //! Construct the Differential Form.
    /*! For Poisson equation, it is confirmed. The Matrix is saved as Sparse Form */
    void makeMatrix();
    //! Construct the RHS.
    /*! Used for Poisson Equation*/
    void makeConstant(double source(double x, double y));
    //! Construct RHS with source vector, time contained.
    /*! Used for Parabola Equation*/
    std::vector<double> makeConstant(std::vector<double> source, double t);
    //! Print the main equation of node (ix, iy)
    void printNodeEQ(int ix, int iy) const;

    //! Operation on PDEMatrix: Get Matrix Element of PDEMatrix
    double GetMatrixElement(int nx, int ny);
    //! Operation on PDEMatrix: Directly change Matrix Element
    void ChangeMatrix(int nx, int ny, double value);
    //! Operation on PDEMatrix: Multiple a factor on the equation of main node (ix, iy).
    /*! This method is used for bad point */
    void ChangeWeight(int ix, int iy, double weight);

    //! Return Nx*Ny
    int GetMatrixSize() const { return Nx * Ny; }
    //! Return PDEMatrix
    const SparseMatrix& GetPDEMatrix() const { return PDEMatrix; }
    //! Reset PDEMatrix
    void SetPDEMatrix(const SparseMatrix& mat) { PDEMatrix = mat; }
    //! Set RHS of PDE
    void setPDEConstant(const std::vector<double>& rhs) { PDEConstant = rhs; }

    //! Form of the Solution of PDE.
    /*! Using Gauss-Sider Iteration Method to solve the PDE problem. */
    SparseMatrix::Solution_info SolvePDE(double tol = 1e-5, int maxIteration = 10000);
};

//! A Boundary Function for Testing.
double TestBoundaryFunc(double x, double y);
//! Write std::vector<double> to a Binary Form file .bin.
void write_vector_to_file(const std::string& filename, const std::vector<double>& vec);




#endif //WEEK1_POISSONEQ_H