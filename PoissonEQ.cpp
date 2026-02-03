//
// Created by 23923 on 2025/10/5.
//

#include "PoissonEQ.h"
#include "sparse_matrix.h"

PolyPoissonEQ::PolyPoissonEQ(int Nx, int Ny, double hx, double hy)
    : Nx(Nx), Ny(Ny), hx(hx), hy(hy), PDEMatrix(SparseMatrix(Nx*Ny,Ny*Ny)), PDEConstant(Nx*Ny,0.0){
};

// 新增：带多边形区域的构造函数
PolyPoissonEQ::PolyPoissonEQ(int Nx, int Ny, double hx, double hy, const PolygonRegion& region, double Boundary_tol)
    : Nx(Nx), Ny(Ny), hx(hx), hy(hy), Boundary_tol(Boundary_tol),
      PDEMatrix(SparseMatrix(Nx*Ny, Nx*Ny)), PDEConstant(Nx*Ny, 0.0),
      polyRegion(new PolygonRegion(region)), isActiveNode(Nx*Ny, false) {
    // 标记活跃节点（多边形内部）
    for (int i = 0; i < Nx*Ny; ++i) {
        double x = Index2PositionX(i);
        double y = Index2PositionY(i);
        isActiveNode[i] = polyRegion->isInside(x, y);
    }
};

double TestBoundaryFunc(double x, double y) {
    return x*x-y*y;
}

double PolyPoissonEQ::Index2PositionX(int index) {
    double posx = hx*(index/Ny);
    return posx;
}

double PolyPoissonEQ::Index2PositionY(int index) {
    double posy = hy*(index%Ny);
    return posy;
}

int PolyPoissonEQ::Index2NodeX(int index) {
    int posx = index%Ny;
    return posx;
}

int PolyPoissonEQ::Index2NodeY(int index) {
    int posy = index/Ny;
    return posy;
}


void PolyPoissonEQ::makeMatrix() {
    double average_weight=10.0;
    std::vector<double> ValList;
    std::vector<int> ColList;
    std::vector<int> RowList;
    RowList.push_back(0);

    for (int ix = 0; ix < Nx; ++ix) {
        for (int iy = 0; iy < Ny; ++iy) {
            int index = ix * Ny + iy;
            double x = Index2PositionX(index);
            double y = Index2PositionY(index);

            // Outside the polyRegion node is set to 0
            if (!isActiveNode[index]) {
                ValList.push_back(1.0);
                ColList.push_back(index);
                RowList.push_back(RowList.back()+1);
                continue;
            }

            // Check if the node is on the boundary
            PolyBoundaryType btype;
            double bvalue, bnx, bny, balpha, bbeta;
            if (polyRegion->isOnBoundary(x, y, btype, bvalue, bnx, bny,balpha,
                bbeta, Boundary_tol)) {
                if (btype == POLY_DIRICHLET) {
                    // Dirichlet BC
                    ValList.push_back(1.0);
                    ColList.push_back(index);
                    RowList.push_back(RowList.back() + 1);
                } else {
                    // Neumann BC
                    int right = (ix+1) * Ny + iy;
                    int up = ix * Ny + iy-1;
                    int left = (ix+1) * Ny + iy;
                    int down = ix * Ny + iy+1;




                    /*// Inside node: Take average of its neighbors
                    std::vector<int> neighbors = {
                        (ix - 1) * Ny + iy,    // left
                        ix * Ny + (iy - 1),    // up
                        ix * Ny + iy,          // self
                        ix * Ny + (iy + 1),    // down
                        (ix + 1) * Ny + iy     // right
                    };
                    int count = 0;
                    int count1 = 0;
                    for (int i = 0; i < 5; ++i) {
                        int neighbor = neighbors[i];
                        PolyBoundaryType new_type;
                        double newvalue, newnx, newny, newalpha, newbeta;
                        // 检查邻点是否都是边界点
                        if (neighbor >= 0 && neighbor < Nx*Ny && isActiveNode[neighbor]) {
                            count ++;
                            if (polyRegion->isOnBoundary(Index2PositionX(neighbor), Index2PositionY(neighbor), new_type,newvalue,newnx,newny
                            ,newalpha, newbeta, Boundary_tol)){
                                count1++;
                            }
                        }
                    }
                    // 若是则该点等于周围点的平均值
                    if (count == count1) {
                        int count2 = 0;
                        int main_index = 0;
                        for (int i = 0; i < 5; ++i) {
                            int neighbor = neighbors[i];
                            PolyBoundaryType new_type2;
                            double newvalue2, newnx2, newny2, newalpha2, newbeta2;
                            if (neighbor >= 0 && neighbor < Nx*Ny && polyRegion->isOnBoundary(x, y, new_type2,newvalue2,newnx2,newny2
                            ,newalpha2, newbeta2, Boundary_tol)) {
                                ValList.push_back(average_weight);
                                if (i == 2) {
                                    main_index = ValList.size();
                                }
                                ColList.push_back(neighbor);
                                count2++;
                            }
                        }
                        ValList[main_index-1] = -average_weight*count2;
                        RowList.push_back(RowList.back() + count2);
                        continue;
                    }*/







                    double bangle_tol = 1e-5;

                    if (ix + 1 < Nx && isActiveNode[right] && balpha+bbeta*(bny / hy+bnx / hx) && bnx>=-bangle_tol) {
                        if (iy - 1 >=0 && isActiveNode[up] && bny >=-bangle_tol) {
                            ValList.push_back(balpha+bbeta*(bny / hy+bnx / hx));
                            ColList.push_back(index);
                            ValList.push_back(-bbeta*bny / hy);  // 法向导数近似
                            ColList.push_back(up);
                            ValList.push_back(-bbeta*bnx / hx);  // 法向导数近似
                            ColList.push_back(right);
                            RowList.push_back(RowList.back() + 3);
                        }else {
                            if (iy+1 < Ny && isActiveNode[down] && balpha+bbeta*(bnx / hx-bny / hy) && bny <= bangle_tol) {
                                ValList.push_back(balpha+bbeta*(bnx / hx-bny / hy));
                                ColList.push_back(index);
                                ValList.push_back(-bbeta*bnx / hx);  // 法向导数近似
                                ColList.push_back(right);
                                ValList.push_back(bbeta*bny / hy);  // 法向导数近似
                                ColList.push_back(down);
                                RowList.push_back(RowList.back() + 3);
                            }else{
                                // Inside node: Take average of its neighbors
                                std::vector<int> neighbors = {
                                    (ix - 1) * Ny + iy,    // left
                                    ix * Ny + (iy - 1),    // up
                                    ix * Ny + iy,          // self
                                    ix * Ny + (iy + 1),    // down
                                    (ix + 1) * Ny + iy     // right
                                };
                                std::vector<double> coeffs = {1.0/(hx*hx), 1.0/(hy*hy), -2.0/(hy*hy)-2.0/(hx*hx), 1.0/(hy*hy), 1.0/(hx*hx)};
                                int count = 0;
                                int main_index = 0;
                                for (int i = 0; i < 5; ++i) {
                                    int neighbor = neighbors[i];
                                    // 检查邻点是否在多边形内
                                    if (neighbor >= 0 && neighbor < Nx*Ny && isActiveNode[neighbor]) {
                                        ValList.push_back(average_weight);
                                        if (i == 2) {
                                            main_index = ValList.size();
                                        }
                                        ColList.push_back(neighbor);
                                        count++;
                                    }
                                }
                                ValList[main_index-1] = -average_weight*count;
                                RowList.push_back(RowList.back() + count);
                            }
                        }
                    } else {
                        // 边界节点无右邻点时，用左邻点
                        if (ix + 1 < Nx && isActiveNode[up] && balpha+bbeta*(bny / hy-bnx / hx)&& bnx<bangle_tol) {
                            ValList.push_back(-bbeta*bny / hy);  // 法向导数近似
                            ColList.push_back(up);
                            ValList.push_back(bbeta*bnx / hx);  // 法向导数近似
                            ColList.push_back(left);
                            ValList.push_back(balpha+bbeta*(bny / hy-bnx / hx));
                            ColList.push_back(index);
                            RowList.push_back(RowList.back() + 3);
                        }else {
                            if (iy+1 < Ny && isActiveNode[down] && balpha-bbeta *(bnx / hx+bny / hy)&& bnx<bangle_tol) {
                                ValList.push_back(bbeta*bnx / hx);  // 法向导数近似
                                ColList.push_back(left);
                                ValList.push_back(balpha-bbeta *(bnx / hx+bny / hy));
                                ColList.push_back(index);
                                ValList.push_back(bbeta*bny / hy);  // 法向导数近似
                                ColList.push_back(down);
                                RowList.push_back(RowList.back() + 3);
                            }else {
                                // Inside node: Take average of its neighbors
                                std::vector<int> neighbors = {
                                    (ix - 1) * Ny + iy,    // left
                                    ix * Ny + (iy - 1),    // up
                                    ix * Ny + iy,          // self
                                    ix * Ny + (iy + 1),    // down
                                    (ix + 1) * Ny + iy     // right
                                };
                                std::vector<double> coeffs = {1.0/(hx*hx), 1.0/(hy*hy), -2.0/(hy*hy)-2.0/(hx*hx), 1.0/(hy*hy), 1.0/(hx*hx)};
                                int main_index = 0;
                                int count = 0;
                                for (int i = 0; i < 5; ++i) {
                                    int neighbor = neighbors[i];
                                    // 检查邻点是否在多边形内
                                    if (neighbor >= 0 && neighbor < Nx*Ny && isActiveNode[neighbor]) {
                                        ValList.push_back(average_weight);
                                        if (i == 2) {
                                            main_index = ValList.size();
                                        }
                                        ColList.push_back(neighbor);
                                        count++;
                                    }
                                }
                                ValList[main_index-1] = -average_weight*count;
                                RowList.push_back(RowList.back() + count);
                            }
                        }
                    }
                }
            } else {
                // Inside node: Take average of its neighbors
                std::vector<int> neighbors = {
                    (ix - 1) * Ny + iy,    // left
                    ix * Ny + (iy - 1),    // up
                    ix * Ny + iy,          // self
                    ix * Ny + (iy + 1),    // down
                    (ix + 1) * Ny + iy     // right
                };
                std::vector<double> coeffs = {1.0/(hx*hx), 1.0/(hy*hy), -2.0/(hy*hy)-2.0/(hx*hx), 1.0/(hy*hy), 1.0/(hx*hx)};
                int count = 0;
                for (int i = 0; i < 5; ++i) {
                    int neighbor = neighbors[i];
                    // 检查邻点是否在多边形内
                    if (neighbor >= 0 && neighbor < Nx*Ny && isActiveNode[neighbor]) {
                        ValList.push_back(coeffs[i]);
                        ColList.push_back(neighbor);
                        count++;
                    }
                }
                RowList.push_back(RowList.back() + count);
            }
        }
    }
    PDEMatrix.setMatrix(ValList, ColList, RowList);
}

// 重写右端项：处理边界条件和源项
void PolyPoissonEQ::makeConstant(double source(double x, double y)) {
    for (int i = 0; i < Nx*Ny; ++i) {
        if (!isActiveNode[i]) {
            PDEConstant[i] = 0.0;
            continue;
        }
        double x = Index2PositionX(i);
        double y = Index2PositionY(i);
        PolyBoundaryType btype;
        double bvalue, bnx, bny, balpha, bbeta;
        if (polyRegion->isOnBoundary(x, y, btype, bvalue, bnx, bny,balpha,
            bbeta, Boundary_tol)) {
            if (btype == POLY_DIRICHLET) {
                PDEConstant[i] = bvalue;      // Dirichlet BC
            } else {
                PDEConstant[i] = bvalue;      // Neumann BC & Mix BC
            }
        } else {
            PDEConstant[i] = source(x, y);    // Inside Node: use source function
        }
    }
}

std::vector<double> PolyPoissonEQ::makeConstant(std::vector<double> source, double t) {
    for (int i = 0; i < Nx*Ny; ++i) {
        if (!isActiveNode[i]) {
            PDEConstant[i] = 0.0;
            continue;
        }
        double x = Index2PositionX(i);
        double y = Index2PositionY(i);
        PolyBoundaryType btype;
        double bvalue, bnx, bny, balpha, bbeta;
        if (polyRegion->isOnBoundary(x, y, t, btype, bvalue, bnx, bny,balpha,
            bbeta, Boundary_tol)) {
            if (btype == POLY_DIRICHLET) {
                PDEConstant[i] = bvalue;      // Dirichlet BC
            } else {
                PDEConstant[i] = bvalue;      // Neumann BC & Mix BC
            }
            } else {
                PDEConstant[i] = source[i];    // Inside Node: use source function
            }
    }
    return PDEConstant;
}


void PolyPoissonEQ::printNodeEQ(int ix, int iy) const {
    PDEMatrix.printRow(ix * Ny + iy);
}

double PolyPoissonEQ::GetMatrixElement(int nx, int ny){
    return PDEMatrix.getValue(nx, ny);
}

void PolyPoissonEQ::ChangeMatrix(int nx, int ny, double value) {
    PDEMatrix.setValue(nx, ny, value);
}

void PolyPoissonEQ::ChangeWeight(int ix, int iy, double weight) {
    PDEMatrix.rowMultiply(ix * Ny + iy, weight);
}


SparseMatrix::Solution_info PolyPoissonEQ::SolvePDE(double tol, int maxIteration){
    return PDEMatrix.gaussSeidel(PDEConstant, tol, maxIteration);
}

void write_vector_to_file(const std::string& filename, const std::vector<double>& vec) {
    // 1. 打开文件流，并设置为二进制写入模式
    std::ofstream outFile(filename, std::ios::binary);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }

    // 2. 首先，将 vector 的大小（元素数量）写入文件
    //    这样 Python 读取时就知道要读取多少个 double。
    size_t size = vec.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // 3. 然后，将 vector 的实际数据写入文件
    //    vec.data() 返回指向底层数据数组的指针 (C++11+)
    //    sizeof(double) * size 是数据的总字节数
    outFile.write(reinterpret_cast<const char*>(vec.data()), sizeof(double) * size);

    std::cout << "Successfully wrote " << vec.size() << " doubles to " << filename << std::endl;

    // 文件流会在 outFile 对象销毁时自动关闭
}
