//
// Created by 23923 on 2025/10/13.
//

#include "ParabolaPDE.h"


ParabolaCN::ParabolaCN(const PolyPoissonEQ& ellipse_solver,
               double alpha, double dt, int total_steps)
        : ellipse_solver(ellipse_solver),
          alpha(alpha), dt(dt), total_steps(total_steps), PDEMatrix(SparseMatrix(ellipse_solver.getNx(),
              ellipse_solver.getNy())) {
    Nx = ellipse_solver.getNx();
    Ny = ellipse_solver.getNy();
    hx = ellipse_solver.gethx();
    hy = ellipse_solver.gethy();
}

std::vector<double> ParabolaCN::setInitialCondition(std::function<double(double, double)> initFunc) {
    std::vector<double> u0(Nx * Ny, 0.0);
    auto  activeNode = ellipse_solver.getActiveNode();
    for (int i = 0; i < Nx * Ny; ++i) {
        if (!activeNode[i]) continue;
        double x = ellipse_solver.Index2PositionX(i);
        double y = ellipse_solver.Index2PositionY(i);
        u0[i] = initFunc(x, y);
    }
    return u0;
}

std::vector<double> ParabolaCN::makeRHS(const std::vector<double>& u_prev,
                                  std::function<double(double, double, double)> sourceFunc,
                                  double t_prev, double t_curr) {
    auto A = ellipse_solver.GetPDEMatrix();
    std::vector<double> sourceTerm;
    for (int i = 0; i < Nx*Ny; i++) {
        double x = ellipse_solver.Index2PositionX(i);
        double y = ellipse_solver.Index2PositionY(i);

        sourceTerm.push_back(dt * sourceFunc(x, y, (t_prev+t_curr)/2));
    }
    vectorAdd(sourceTerm,1.0,u_prev);
    vectorAdd(sourceTerm,-alpha*dt/2,A.multiplyVector(u_prev));
    return  ellipse_solver.makeConstant(sourceTerm, (t_prev+t_curr)/2);
}

void ParabolaCN::makePDE() {
    ellipse_solver.makeMatrix();
    PDEMatrix = ellipse_solver.GetPDEMatrix();
    std::vector<int> rowList = PDEMatrix.get_row_ptr();
    std::vector<int> colList = PDEMatrix.get_col_indices();

    for (int i = 0; i < Nx*Ny; i++) {
        if (rowList[i+1]-rowList[i] > 3) {
            PDEMatrix.rowMultiply(i, alpha*dt/2);
            double diag = PDEMatrix.getValue(i, i);
            PDEMatrix.setValue(i, i, diag+1.0);
        }
        if (PDEMatrix.getValue(i, i)==0) {
            PDEMatrix.printRow(i);
        }
    }
}

void ParabolaCN::solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc, std::string filename) {
    makePDE();
    // Initialize
    std::vector<double> u_prev = setInitialCondition(initFunc);
    saveStepResult(u_prev, 0,filename);

    std::cout << "CN form solver started, totally:" << total_steps << " steps" << std::endl;

    // Iteration
    for (int step = 0; step < total_steps; ++step) {
        double t_prev = step * dt;       // 上一时刻
        double t_curr = (step + 1) * dt; // 当前时刻

        // 1. Make RHS
        std::vector<double> rhs = makeRHS(u_prev, sourceFunc, t_prev, t_curr);

        // 2. Gauss-Sider Solver
        auto solution_t = PDEMatrix.gaussSeidel(rhs,1e-5, 5000);
        std::vector<double> u_curr = solution_t.solution;  // 获取当前时刻解

        // 3. Update result and save.
        u_prev = u_curr;
        saveStepResult(u_curr, step + 1, filename);
        // Output State
        if ((step + 1) % 10 == 0) {
            std::cout << "Finish " << step + 1 << "/" << total_steps
                      << ". Residual:" << solution_t.residual << std::endl;
        }
    }

    std::cout << "CN finished!" << std::endl;
}

std::vector<double> ParabolaCN::solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc) {
    makePDE();
    // Initialize
    std::vector<double> u_prev = setInitialCondition(initFunc);

    std::cout << "CN form solver started, totally:" << total_steps << " steps" << std::endl;

    // Iteration
    for (int step = 0; step < total_steps; ++step) {
        double t_prev = step * dt;       // 上一时刻
        double t_curr = (step + 1) * dt; // 当前时刻

        // 1. Make RHS
        std::vector<double> rhs = makeRHS(u_prev, sourceFunc, t_prev, t_curr);

        // 2. Gauss-Sider Solver
        auto solution_t = PDEMatrix.gaussSeidel(rhs,1e-5, 5000);
        std::vector<double> u_curr = solution_t.solution;  // 获取当前时刻解

        // 3. Update result and save.
        u_prev = u_curr;
        // Output State
        if ((step + 1) % 10 == 0) {
            std::cout << "Finish " << step + 1 << "/" << total_steps
                      << ". Residual:" << solution_t.residual << std::endl;
        }
    }

    std::cout << "CN finished!" << std::endl;
    return u_prev;
}

void ParabolaCN::saveStepResult(const std::vector<double>& u, int step, std::string filename) const {
    filename = filename + std::to_string(step) + ".bin";
    write_vector_to_file(filename, u);  // 复用椭圆求解器的二进制写入函数
}


ParabolaExplicit::ParabolaExplicit(const PolyPoissonEQ& ellipse_solver,
               double alpha, double dt, int total_steps)
        : ellipse_solver(ellipse_solver),
          alpha(alpha), dt(dt), total_steps(total_steps), PDEMatrix(SparseMatrix(ellipse_solver.getNx(),
              ellipse_solver.getNy())) {
    Nx = ellipse_solver.getNx();
    Ny = ellipse_solver.getNy();
    hx = ellipse_solver.gethx();
    hy = ellipse_solver.gethy();
}

std::vector<double> ParabolaExplicit::setInitialCondition(std::function<double(double, double)> initFunc) {
    std::vector<double> u0(Nx * Ny, 0.0);
    auto  activeNode = ellipse_solver.getActiveNode();
    for (int i = 0; i < Nx * Ny; ++i) {
        if (!activeNode[i]) continue;
        double x = ellipse_solver.Index2PositionX(i);
        double y = ellipse_solver.Index2PositionY(i);
        u0[i] = initFunc(x, y);
    }
    return u0;
}

std::vector<double> ParabolaExplicit::makeRHS(const std::vector<double>& u_prev,
                                  std::function<double(double, double, double)> sourceFunc,
                                  double t_prev) {
    auto A = ellipse_solver.GetPDEMatrix();
    std::vector<double> sourceTerm;
    for (int i = 0; i < Nx*Ny; i++) {
        double x = ellipse_solver.Index2PositionX(i);
        double y = ellipse_solver.Index2PositionY(i);

        sourceTerm.push_back(dt * sourceFunc(x, y, t_prev));
    }
    vectorAdd(sourceTerm,1.0,u_prev);
    vectorAdd(sourceTerm,-alpha*dt,A.multiplyVector(u_prev));
    return  ellipse_solver.makeConstant(sourceTerm, t_prev);
}

void ParabolaExplicit::makePDE() {
    ellipse_solver.makeMatrix();
    PDEMatrix = ellipse_solver.GetPDEMatrix();
    std::vector<int> rowList = PDEMatrix.get_row_ptr();
    std::vector<int> colList = PDEMatrix.get_col_indices();

    for (int i = 0; i < Nx*Ny; i++) {
        if (rowList[i+1]-rowList[i] > 3) {
            PDEMatrix.rowMultiply(i, 0);
            double diag = PDEMatrix.getValue(i, i);
            PDEMatrix.setValue(i, i, 1.0);
        }
        if (PDEMatrix.getValue(i, i)==0) {
            PDEMatrix.printRow(i);
        }
    }
}

void ParabolaExplicit::solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc, std::string filename) {
    makePDE();
    // Initialize
    std::vector<double> u_prev = setInitialCondition(initFunc);
    saveStepResult(u_prev, 0,filename);

    std::cout << "CN form solver started, totally:" << total_steps << " steps" << std::endl;

    // Iteration
    for (int step = 0; step < total_steps; ++step) {
        double t_prev = step * dt;       // previous time

        // 1. Make RHS
        std::vector<double> rhs = makeRHS(u_prev, sourceFunc, t_prev);

        // 2. Gauss-Sider Solver
        auto solution_t = PDEMatrix.gaussSeidel(rhs,1e-5, 5000);
        std::vector<double> u_curr = solution_t.solution;  // 获取当前时刻解

        // 3. Update result and save.
        u_prev = u_curr;
        saveStepResult(u_curr, step + 1, filename);
        // Output State
        if ((step + 1) % 10 == 0) {
            std::cout << "Finish " << step + 1 << "/" << total_steps
                      << ". Residual:" << solution_t.residual << std::endl;
        }
    }

    std::cout << "CN finished!" << std::endl;
}

std::vector<double> ParabolaExplicit::solve(std::function<double(double, double)> initFunc,
               std::function<double(double, double, double)> sourceFunc) {
    makePDE();
    // Initialize
    std::vector<double> u_prev = setInitialCondition(initFunc);

    std::cout << "CN form solver started, totally:" << total_steps << " steps" << std::endl;

    // Iteration
    for (int step = 0; step < total_steps; ++step) {
        double t_prev = step * dt;       // Previous time

        // 1. Make RHS
        std::vector<double> rhs = makeRHS(u_prev, sourceFunc, t_prev);

        // 2. Gauss-Sider Solver
        auto solution_t = PDEMatrix.gaussSeidel(rhs,1e-4, 5000);
        std::vector<double> u_curr = solution_t.solution;  // 获取当前时刻解

        // 3. Update result and save.
        u_prev = u_curr;
        // Output State
        if ((step + 1) % 10 == 0) {
            std::cout << "Finish " << step + 1 << "/" << total_steps
                      << ". Residual:" << solution_t.residual << std::endl;
        }
    }

    std::cout << "CN finished!" << std::endl;
    return u_prev;
}

void ParabolaExplicit::saveStepResult(const std::vector<double>& u, int step, std::string filename) const {
    filename = filename + std::to_string(step) + ".bin";
    write_vector_to_file(filename, u);  // 复用椭圆求解器的二进制写入函数
}