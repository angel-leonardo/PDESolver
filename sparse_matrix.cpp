#include "sparse_matrix.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <utility>





// 构造函数
SparseMatrix::SparseMatrix(int rows, int cols)
    : rows(rows), cols(cols), non_zero_count(0) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("Rows and cols must be positive");
    }
    row_ptr.resize(rows + 1, 0);
}

// 拷贝构造函数
SparseMatrix::SparseMatrix(const SparseMatrix& other)
    : rows(other.rows), cols(other.cols), non_zero_count(other.non_zero_count),
      values(other.values), col_indices(other.col_indices), row_ptr(other.row_ptr) {}

// 赋值运算符重载
SparseMatrix& SparseMatrix::operator=(const SparseMatrix& other) {
    if (this != &other) {
        rows = other.rows;
        cols = other.cols;
        non_zero_count = other.non_zero_count;
        values = other.values;
        col_indices = other.col_indices;
        row_ptr = other.row_ptr;
    }
    return *this;
}

// 设置元素值
void SparseMatrix::setValue(int row, int col, double value) {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Element index out of range");
    }

    // 如果值为0，不需要存储
    if (value == 0.0) {
        // 这里可以添加删除现有非零元素的逻辑
        return;
    }

    // 查找行中是否已存在该列的元素
    int start = row_ptr[row];
    int end = row_ptr[row + 1];
    int pos = end;

    for (int i = start; i < end; ++i) {
        if (col_indices[i] < col) {
            continue;
        }
        if (col_indices[i] == col) {
            // 更新已有元素
            values[i] = value;
            return;
        }
        pos=i;
        break;
    }

    // 插入新元素
    values.insert(values.begin() + pos, value);
    col_indices.insert(col_indices.begin() + pos, col);

    // 更新行指针
    for (int i = row + 1; i <= rows; ++i) {
        row_ptr[i]++;
    }

    non_zero_count++;
}

// 直接设置压缩矩阵格式
void SparseMatrix::setMatrix(std::vector<double> value_l, std::vector<int> col_l, std::vector<int> row_l) {
    values = std::move(value_l);
    col_indices = std::move(col_l);
    row_ptr = std::move(row_l);
    non_zero_count = row_ptr[rows];
}

// 获取元素值
double SparseMatrix::getValue(int row, int col) const {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Element index out of range");
    }

    // 在该行的非零元素中查找
    int start = row_ptr[row];
    int end = row_ptr[row + 1];

    for (int i = start; i < end; ++i) {
        if (col_indices[i] == col) {
            return values[i];
        }
    }

    // 未找到，返回0
    return 0.0;
}

// 矩阵转置
SparseMatrix SparseMatrix::transpose() const {
    SparseMatrix result(cols, rows);

    // 统计每列的非零元素数量
    std::vector<int> col_counts(cols, 0);
    for (int i = 0; i < non_zero_count; ++i) {
        col_counts[col_indices[i]]++;
    }

    // 计算转置矩阵的行指针
    result.row_ptr[0] = 0;
    for (int i = 0; i < cols; ++i) {
        result.row_ptr[i + 1] = result.row_ptr[i] + col_counts[i];
    }

    // 重置计数，用于定位插入位置
    std::fill(col_counts.begin(), col_counts.end(), 0);
    result.values.resize(non_zero_count, 0.0);
    result.col_indices.resize(non_zero_count, 0);
    // 填充转置矩阵
    for (int i = 0; i < rows; ++i) {
        int start = row_ptr[i];
        int end = row_ptr[i + 1];

        for (int j = start; j < end; ++j) {
            int col = col_indices[j];
            double val = values[j];

            // 计算在转置矩阵中的位置
            int pos = result.row_ptr[col] + col_counts[col];
            col_counts[col]++;


            result.values[pos] = val;
            result.col_indices[pos] = i;
        }
    }

    result.non_zero_count = non_zero_count;
    return result;
}

//! Matrix add
SparseMatrix SparseMatrix::add(const SparseMatrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimension mismatch");
    }

    SparseMatrix result(rows, cols);
    // 预分配空间（可选优化，减少动态扩容开销）
    result.values.reserve(non_zero_count + other.non_zero_count);
    result.col_indices.reserve(non_zero_count + other.non_zero_count);
    result.row_ptr.reserve(rows + 1);
    result.row_ptr.push_back(0);  // 首行指针固定为0

    for (int i = 0; i < rows; ++i) {
        // 当前行A的非零元素范围
        int startA = row_ptr[i];
        int endA = row_ptr[i + 1];
        // 当前行B的非零元素范围
        int startB = other.row_ptr[i];
        int endB = other.row_ptr[i + 1];

        int pA = startA;  // A的当前元素指针
        int pB = startB;  // B的当前元素指针

        // 合并当前行的非零元素
        while (pA < endA && pB < endB) {
            int colA = col_indices[pA];
            int colB = other.col_indices[pB];

            if (colA == colB) {
                // 同列元素，相加（非零才保留）
                double sum_val = values[pA] + other.values[pB];
                result.setValue(i, colA, sum_val);

                pA++;
                pB++;

            } else if (colA > colB) {
                // B有此列元素，A无
                result.setValue(i, colB, other.values[pB]);
                pB++;
            } else {
                // A有此列元素，B无
                result.setValue(i, colA, values[pA]);
                pA++;
            }
        }

        // 处理A中剩余的非零元素
        while (pA < endA) {
            result.setValue(i, col_indices[pA], values[pA]);
            pA++;
        }

        // 处理B中剩余的非零元素
        while (pB < endB) {
            result.setValue(i, other.col_indices[pB], other.values[pB]);
            pB++;
        }

        // 更新当前行的结束指针
        result.row_ptr.push_back(result.values.size());
    }

    // 设置非零元素总数
    result.non_zero_count = result.values.size();
    return result;
}

//! Matrix Multiply
SparseMatrix SparseMatrix::multiply(const SparseMatrix& other) const {
    if (cols != other.rows) {
        throw std::invalid_argument("Matrix dimension mismatch");
    }

    SparseMatrix result(rows, other.cols);
    SparseMatrix otherTransposed = other.transpose();  // 转置以提高缓存利用率

    // 计算每个元素
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            double sum = 0.0;

            // 获取当前行的非零元素范围
            int start1 = row_ptr[i];
            int end1 = row_ptr[i + 1];

            // 获取转置矩阵中当前行（即原矩阵的列）的非零元素范围
            int start2 = otherTransposed.row_ptr[j];
            int end2 = otherTransposed.row_ptr[j + 1];

            // 计算点积
            int p1 = start1, p2 = start2;
            while (p1 < end1 && p2 < end2) {
                int col1 = col_indices[p1];
                int row2 = otherTransposed.col_indices[p2];  // 原矩阵的行索引

                if (col1 == row2) {
                    sum += values[p1] * otherTransposed.values[p2];
                    p1++;
                    p2++;
                } else if (col1 < row2) {
                    p1++;
                } else {
                    p2++;
                }
            }

            if (sum != 0.0) {
                result.setValue(i, j, sum);
            }
        }
    }

    return result;
}

//! Print Matrix in matrix form
void SparseMatrix::printDense() const {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << getValue(i, j) << "\t";
        }
        std::cout << std::endl;
    }
}
//! Print Matrix in RCS
void SparseMatrix::printSparse() const {
    std::cout << "row: " << rows << ", col: " << cols << ", elements: " << non_zero_count << std::endl;

    std::cout << "values: ";
    for (double v : values) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::cout << "col_indices: ";
    for (int c : col_indices) {
        std::cout << c << " ";
    }
    std::cout << std::endl;

    std::cout << "row_ptr: ";
    for (int r : row_ptr) {
        std::cout << r << " ";
    }
    std::cout << std::endl;
}
//! Print a row
void SparseMatrix::printRow(int row) const {
    std::cout << "row: " << row << std::endl;
    std::cout << "values: ";
    for (int i= row_ptr[row]; i < row_ptr[row + 1]; ++i) {
        double v = values[i];
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::cout << "col_indices: ";
    for (int i= row_ptr[row]; i < row_ptr[row + 1]; ++i) {
        int c = col_indices[i];
        std::cout << c << " ";
    }
    std::cout << std::endl;
}

std::vector<double> SparseMatrix::multiplyVector(const std::vector<double>& vec) const {
    if (vec.size() != cols) {
        throw std::invalid_argument("Vector size does not match matrix columns");
    }

    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        int start = row_ptr[i];
        int end = row_ptr[i + 1];
        for (int j = start; j < end; ++j) {
            result[i] += values[j] * vec[col_indices[j]];
        }
    }
    return result;
}


//! Calculate ||vector||2
double norm(const std::vector<double>& vec) {
    double sum = 0.0;
    for (double v : vec) {
        sum += v * v;
    }
    return std::sqrt(sum);
}

//! Calculate vector minus
void vectorSubtract(std::vector<double>& a, const std::vector<double>& b) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] -= b[i];
    }
}

//! Calculate vector plus
void vectorAdd(std::vector<double>& a, double alpha, const std::vector<double>& b) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] += alpha * b[i];
    }
}

//! Calculate vector number multiply
void vectorNumMultiply(std::vector<double>& a, double alpha) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = alpha * a[i];
    }
}

//! Calculate vector inner product
double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}


//! Elementary Transformation
//! rowMultiply
void SparseMatrix::rowMultiply(int row, double factor){
    int i = row_ptr[row];
    int end = row_ptr[row + 1];
    while (i < end) {
        values[i] = values[i] * factor;
        i++;
    }
}
//! rowMultiAdd
void SparseMatrix::rowMultiAdd(int row1, double factor, int row2) {
    for (int i = row_ptr[row1]; i < row_ptr[row1+1]; ++i) {
        double temp = getValue(row2, col_indices[i]);
        setValue(row2, col_indices[i], temp + getValue(row1, col_indices[i])*factor);
    }
};


// Gauss-Seidel迭代法
SparseMatrix::Solution_info SparseMatrix::gaussSeidel(const std::vector<double>& b,
                                                      double tol,
                                                      int maxIter) const {

    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square for Gauss-Seidel iteration");
    }
    if (b.size() != rows) {
        throw std::invalid_argument("Vector size does not match matrix dimensions");
    }

    std::vector<double> x(rows, 0.0);  // 初始解向量
    std::vector<double> x_old(rows);
    Solution_info solution;
    solution.iteration = 0;
    solution.residual = 100.0;
    do {
        x_old = x;

        // 迭代更新每个分量
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            double diag = 0.0;

            // 计算该行的总和，区分对角线元素
            int start = row_ptr[i];
            int end = row_ptr[i + 1];
            for (int j = start; j < end; ++j) {
                int col = col_indices[j];
                if (col == i) {
                    diag = values[j];  // 对角线元素
                } else {
                    sum += values[j] * x[col];  // 使用最新的x值
                }
            }

            if (std::fabs(diag) < 1e-12) {
                throw std::runtime_error("Matrix has zero diagonal element, Gauss-Seidel cannot proceed");
            }

            // 更新x[i]
            x[i] = (b[i] - sum) / diag;
        }

        // 计算残差
        std::vector<double> Ax = multiplyVector(x);
        vectorSubtract(Ax, b);
        solution.residual = norm(Ax);

        if (solution.residual < tol) {
            break;
        }

    } while (++solution.iteration < maxIter);

    if (solution.iteration >= maxIter) {
        std::cerr << "Gauss-Seidel did not converge within " << maxIter << " iterations. Residual" << solution.residual<< std::endl;
    }

    solution.solution = x;


    return solution;
}

// 共轭梯度法（适用于对称正定矩阵）
// 修正后的共轭梯度法（适用于对称正定矩阵）
SparseMatrix::Solution_info SparseMatrix::conjugateGradient(const std::vector<double>& b,
                                                  double tol,
                                                  int maxIter) const {
    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square for conjugate gradient method");
    }
    if (b.size() != rows) {
        throw std::invalid_argument("Vector size does not match matrix dimensions");
    }
    Solution_info solution;
    solution.iteration = 0;
    solution.residual = 100.0;
    std::vector<double> x(rows, 0.0);  // 初始解向量
    std::vector<double> r = b;         // 初始残差 r0 = b - Ax0 (x0=0)
    std::vector<double> p = r;         // 初始搜索方向
    double r_dot_r = dotProduct(r, r); // r0·r0

    // 检查初始残差
    solution.residual = std::sqrt(r_dot_r);
    if (solution.residual < tol) {
        solution.solution = x;
        return solution;
    }

    for (solution.iteration = 0; solution.iteration < maxIter; ++solution.iteration) {
        std::vector<double> Ap = multiplyVector(p);  // 计算Ap
        double p_dot_Ap = dotProduct(p, Ap);

        // 检查矩阵正定性
        if (p_dot_Ap <= 0.0) {
            throw std::runtime_error("Matrix is not positive definite, conjugate gradient cannot proceed");
        }

        double alpha = r_dot_r / p_dot_Ap;  // 计算步长

        // 更新解 x = x + alpha*p
        vectorAdd(x, alpha, p);

        // 保存旧残差用于计算beta
        std::vector<double> r_old = r;

        // 更新残差 r = r - alpha*Ap
        vectorAdd(r, -alpha, Ap);

        // 计算新残差的点积
        double r_new_dot_r_new = dotProduct(r, r);
        solution.residual = std::sqrt(r_new_dot_r_new);

        // 检查收敛
        if (solution.residual < tol) {

            break;
        }

        // 计算新的搜索方向
        double beta = r_new_dot_r_new / r_dot_r;
        // p_new = r + beta*p_old
        std::vector<double> p_new = r;
        vectorAdd(p_new, beta, p);
        p = p_new;

        r_dot_r = r_new_dot_r_new;
    }

    if (solution.iteration >= maxIter) {
        double final_residual = std::sqrt(dotProduct(r, r));
        std::cerr << "Conjugate gradient did not converge within " << maxIter
                  << " iterations. Final residual: " << final_residual << std::endl;
    }
    solution.solution = x;
    return solution;


}

// 获取对角线部分（i == j）
SparseMatrix SparseMatrix::diagonal() const {
    SparseMatrix result(rows, cols);  // 结果矩阵维度与原矩阵相同

    // 遍历每一行的非零元素
    for (int i = 0; i < rows; ++i) {
        int start = row_ptr[i];
        int end = row_ptr[i + 1];

        for (int j = start; j < end; ++j) {
            int col = col_indices[j];
            double val = values[j];

            // 只保留行索引等于列索引的元素
            if (i == col) {
                result.setValue(i, col, val);
            }
        }
    }

    return result;
}

// 获取上三角部分（i <= j，包含对角线）
SparseMatrix SparseMatrix::upperTriangle() const {
    SparseMatrix result(rows, cols);

    for (int i = 0; i < rows; ++i) {
        int start = row_ptr[i];
        int end = row_ptr[i + 1];

        for (int j = start; j < end; ++j) {
            int col = col_indices[j];
            double val = values[j];

            // 保留行索引小于等于列索引的元素
            if (i < col) {
                result.setValue(i, col, val);
            }
        }
    }

    return result;
}

// 获取下三角部分（i >= j，包含对角线）
SparseMatrix SparseMatrix::lowerTriangle() const {
    SparseMatrix result(rows, cols);

    for (int i = 0; i < rows; ++i) {
        int start = row_ptr[i];
        int end = row_ptr[i + 1];

        for (int j = start; j < end; ++j) {
            int col = col_indices[j];
            double val = values[j];

            // 保留行索引大于等于列索引的元素
            if (i > col) {
                result.setValue(i, col, val);
            }
        }
    }

    return result;
}

SparseMatrix Identity(int dim) {
    SparseMatrix result(dim, dim);
    for (int i = 0; i < dim; ++i) {
        result.setValue(i, i, 1.0);
    }
    return result;
}

double sum(std::vector<double> a) {
    double result = 0.0;
    for (int i = 0; i < a.size(); ++i) {
        result += a[i];
    }
    return result;
}

std::vector<double> vectorIdentity(int num, double value) {
    std::vector<double> result;
    for (int i = 0; i < num; ++i) {
        result.push_back(value);
    }
    return result;
}
