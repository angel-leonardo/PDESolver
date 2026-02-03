//
// Created by 23923 on 2026/1/11.
//

#ifndef VISUALIZEMAIN_PY_HYPERBOLICEQ_H
#define VISUALIZEMAIN_PY_HYPERBOLICEQ_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdint> // 用于size_t（对应Python的struct 'Q'）

// 保存vector<double>到二进制文件（与Python read_vector_from_file兼容）
void save_vector_to_file(const std::vector<double>& vec, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }
    // 写入vector大小（size_t类型，对应Python的struct.unpack('Q')）
    size_t size = vec.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
    // 写入数据
    outFile.write(reinterpret_cast<const char*>(vec.data()), sizeof(double) * size);
    outFile.close();
}

// Lax-Wendroff格式求解一维线性对流方程：∂u/∂t + a*∂u/∂x = 0
void solve_hyperbolic_lax_wendroff() {
    // ===================== 网格/参数配置 =====================
    const double x_start = 0.0;    // x轴起始
    const double x_end = 2.0;      // x轴终止
    const int nx = 201;            // x方向节点数
    const double dx = (x_end - x_start) / (nx - 1); // x步长

    const double t_total = 1.0;    // 总模拟时间
    const int nt = 500;            // 时间步数
    const double dt = t_total / nt; // 时间步长

    const double a = 1.0;          // 对流速度（常数）
    const double cfl = a * dt / dx;// CFL数（稳定性条件：|cfl| ≤ 1）

    // 稳定性检查
    if (fabs(cfl) > 1.0) {
        std::cerr << "警告：CFL数 = " << cfl << " > 1，Lax-Wendroff格式不稳定！" << std::endl;
        // 可选：自动调整dt以满足稳定性
        // dt = dx * 1.0 / fabs(a);
        // nt = static_cast<int>(t_total / dt);
        return;
    }

    // ===================== 初始条件 =====================
    // 示例：高斯脉冲初始条件 u(x,0) = exp(-50*(x-0.5)^2)
    std::vector<double> u(nx, 0.0);
    for (int j = 0; j < nx; ++j) {
        double x = x_start + j * dx;
        u[j] = exp(-50.0 * pow(x - 0.5, 2));
    }

    // 保存初始时间步（t=0）
    save_vector_to_file(u, "hyperbolic_lw_step_0.bin");

    // ===================== Lax-Wendroff时间步进 =====================
    std::vector<double> u_new(nx, 0.0);
    for (int n = 1; n <= nt; ++n) { // 时间步从1到nt
        // 内部节点（j=1到nx-2）
        for (int j = 1; j < nx - 1; ++j) {
            // Lax-Wendroff格式核心公式
            u_new[j] = u[j]
                      - 0.5 * cfl * (u[j+1] - u[j-1])
                      + 0.5 * cfl * cfl * (u[j+1] - 2*u[j] + u[j-1]);
        }

        // 边界条件：周期性边界（对流方程常用）
        u_new[0] = u_new[nx-2];    // 左边界 = 倒数第二个节点
        u_new[nx-1] = u_new[1];    // 右边界 = 第二个节点

        // 更新u为新时间步的值
        u.swap(u_new);

        // 保存当前时间步结果
        std::string filename = "hyperbolic_lw_step_" + std::to_string(n) + ".bin";
        save_vector_to_file(u, filename);

        // 打印进度
        if (n % 50 == 0) {
            std::cout << "已完成时间步: " << n << "/" << nt << std::endl;
        }
    }

    std::cout << "双曲型方程求解完成！共生成 " << nt + 1 << " 个时间步文件。" << std::endl;
}


#endif //VISUALIZEMAIN_PY_HYPERBOLICEQ_H