import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import struct
import os


# 读取C++输出的二进制向量文件
def read_vector_from_file(filename):
    """读取C++写入的std::vector<double>二进制文件"""
    try:
        with open(filename, 'rb') as inFile:
            # 读取向量大小（size_t类型）
            size_bytes = inFile.read(struct.calcsize('Q'))
            if not size_bytes:
                return None
            size = struct.unpack('Q', size_bytes)[0]

            # 读取数据
            data_bytes = inFile.read(struct.calcsize('d') * size)
            if len(data_bytes) != struct.calcsize('d') * size:
                print(f"警告: 文件 {filename} 损坏或不完整")
                return None

            return np.array(struct.unpack('d' * size, data_bytes))

    except FileNotFoundError:
        print(f"错误: 找不到文件 {filename}")
        return None
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return None


# 创建抛物型方程结果的动画可视化
def create_parabola_animation():
    # 网格参数（与C++代码保持一致）
    nx, ny = 101, 201  # 节点数量
    hx, hy = 0.02, 0.02  # 网格步长
    total_steps = 1001  # 总时间步数

    # 读取所有时间步的结果数据
    all_data = []
    for i in range(total_steps):
        filename = f"parabola_cn_step_{i}.bin"
        data = read_vector_from_file(filename)
        if data is None:
            print(f"无法读取文件 {filename}，程序退出")
            return
        # 转换为2D矩阵
        all_data.append(data.reshape((nx, ny)))

    # 确定整个序列的颜色范围（确保颜色映射一致）
    all_values = np.concatenate([data.flatten() for data in all_data])
    vmin, vmax = all_values.min(), all_values.max()
    norm = Normalize(vmin=vmin, vmax=vmax)

    # 创建坐标网格
    x = np.linspace(0, (nx - 1) * hx, nx)
    y = np.linspace(0, (ny - 1) * hy, ny)

    # 创建画布
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.suptitle('抛物型方程数值解时间演化', fontsize=16)

    # 初始化热图
    im = ax.imshow(all_data[0], cmap='viridis', origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()],
                   aspect='auto', norm=norm)

    # 添加颜色条和标签
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('函数值 U(x,y,t)')

    ax.set_xlabel('X坐标')
    ax.set_ylabel('Y坐标')
    time_text = ax.text(0.05, 0.95, f'时间步: 0/{total_steps - 1}',
                        transform=ax.transAxes,
                        bbox=dict(facecolor='white', alpha=0.8))

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # 动画更新函数
    def update(frame):
        im.set_data(all_data[frame])
        time_text.set_text(f'时间步: {frame}/{total_steps - 1}')
        return im, time_text

    # 创建动画
    anim = animation.FuncAnimation(
        fig, update,
        frames=total_steps,
        interval=100,  # 每帧间隔100ms
        blit=True
    )

    # 保存动画（尝试多种编码器兼容）
    try:
        # 尝试保存为GIF
        writer = animation.PillowWriter(fps=10)
        anim.save('parabola_animation.gif', writer=writer, dpi=100)
        print("动画已保存为 'parabola_animation.gif'")
    except Exception as e:
        print(f"保存GIF时出错: {e}")
        try:
            # 尝试保存为MP4（需要ffmpeg）
            anim.save('parabola_animation.mp4', writer='ffmpeg', fps=10, dpi=100)
            print("动画已保存为 'parabola_animation.mp4'")
        except Exception as e2:
            print(f"保存MP4时出错: {e2}")

    # 显示动画
    plt.show()


if __name__ == "__main__":
    create_parabola_animation()
