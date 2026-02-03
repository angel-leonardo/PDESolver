#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <iostream>


//!  A Sparse Matrix in RCS.
/*!
 * This is a Sparse Matrix in RCS. With these data saved:
 *  Rows of Matrix
 *   Column of Matrix
 *   Nonzero element count
 *   List of value
 *   List of column indices
 *   Start point of index of each row
*/
class SparseMatrix {
    //! Form of Sparse Matrix.
    /*! More detailed enum description. */
private:
    /*! Rows of Matrix */
    int rows;
    /*! Column of Matrix */
    int cols;
    /*! Nonzero element count */
    int non_zero_count;


    /*! List of value */
    std::vector<double> values;
    /*! List of column indices */
    std::vector<int> col_indices;
    /*! Start point of index of each row */
    std::vector<int> row_ptr;

public:
    //! Construct a SparseMatrix.
    /*! Set rows and columns of Sparse Matrix */
    SparseMatrix(int rows, int cols);

    //! Copy of SparseMatrix
    SparseMatrix(const SparseMatrix& other);
    
    // 析构函数
    ~SparseMatrix() = default;
    
    //! Reload "=" of SparseMatrix
    SparseMatrix& operator=(const SparseMatrix& other);
    
    //! Get row information of SparseMatrix
    int getRows() const { return rows; }
    //! Get column information of SparseMatrix
    int getCols() const { return cols; }
    //! Get nonzero information of SparseMatrix
    int getNonZeroCount() const { return non_zero_count; }
    //! Get row_ptr of SparseMatrix
    std::vector<int> get_row_ptr() const { return row_ptr; }
    //! Get col_indices of SparseMatrix
    std::vector<int> get_col_indices() const { return col_indices; }
    
    //! Insert Value
    /*! Caution: complexity equals to the number of elements. Don't use it in large scale case. */
    void setValue(int row, int col, double value);
    //! Set SparseMatrix Value
    /*! If you have the form of the SparseMatrix, this method is recommended. */
    void setMatrix(std::vector<double> value_l, std::vector<int> col_l, std::vector<int> row_l);
    
    //! Get value of a specific row and column.
    double getValue(int row, int col) const;
    
    //! Transpose of a SparseMatrix
    SparseMatrix transpose() const;


    //! Elementary transformation: Multiply a row with a factor.
    void rowMultiply(int row, double factor);
    //! Elementary transformation: Multiply row1 with factor and add the result to row2.
    void rowMultiAdd(int row1, double factor, int row2);


    //! Get diagonal part of the SparseMatrix
    SparseMatrix diagonal() const;
    //! Upper Triangle part of the SparseMatrix
    /*! Diagonal part not included.*/
    SparseMatrix upperTriangle() const;
    //! Lower Triangle part of the SparseMatrix
    /*! Diagonal part not included.*/
    SparseMatrix lowerTriangle() const;


    //! Add of SparseMatrix
    SparseMatrix add(const SparseMatrix& other) const;
    //! Multiply of SparseMatrix
    SparseMatrix multiply(const SparseMatrix& other) const;


    //! Print Matrix in Dense Form
    void printDense() const;
    //! Print Matrix in Sparse Form
    void printSparse() const;
    //! Print a row of SparseMatrix
    void printRow(int row) const;

    //! Multiply of SparseMatrix to Vector
    std::vector<double> multiplyVector(const std::vector<double>& vec) const;


    //! Multiply the SparseMatrix by a scalar
    SparseMatrix scalarMultiply(double factor) const {
        SparseMatrix result(rows, cols);
        result.values = values;
        for (auto& v : result.values) v *= factor;
        result.col_indices = col_indices;
        result.row_ptr = row_ptr;
        result.non_zero_count = non_zero_count;
        return result;
    }

    //! Form of solution.
    /*! (int iteration, double residual, std::vector<double> solution)
     *  All solution of SparseMatrix is returned in this form.
     *  Iteration will be set as default 0 if it is not a iterative algorithm.
     *  If the algorithm does not converge, will throw an error, and return
     *  residual nan and solution nan.
     */
    struct Solution_info {
        //! Iteration of Solving
        int iteration;
        //! Residual of Solution
        double residual;
        //! Solution vector
        std::vector<double> solution;
    };


    //! Gauss-Seidel iterative algorithm to solve linear equation.
     Solution_info gaussSeidel(const std::vector<double>& b,
                                   double tol = 1e-6,
                                   int maxIter = 1000) const;

    //! Conjugate gradient algorithm to solve linear equation
    Solution_info conjugateGradient(const std::vector<double>& b,
                                         double tol = 1e-5,
                                         int maxIter = 1000) const;
};

//! Calculate ||vector||2
double norm(const std::vector<double>& vec);
//! Calculate vector minus
void vectorSubtract(std::vector<double>& a, const std::vector<double>& b);
//! Calculate vector plus
void vectorAdd(std::vector<double>& a, double alpha, const std::vector<double>& b);
//! Calculate vector number multiply
void vectorNumMultiply(std::vector<double>& a, double alpha);
//! Calculate vector inner product
double dotProduct(const std::vector<double>& a, const std::vector<double>& b);
//! Initialize a SparseMatrix as Identity
SparseMatrix Identity(int dim);
double sum(std::vector<double> a);
std::vector<double> vectorIdentity(int num, double value);

#endif // SPARSE_MATRIX_H
