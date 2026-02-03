import matplotlib.pyplot as plt


def hotimage(two_d_matrix, size= (8,6)):
    """
    Plot Matrix as Hot Image
    :param two_d_matrix: Data of Hot Image
    :param size: Size of Image
    :return: void
    """
    plt.figure(figsize=size)

    # 使用 imshow 绘制热图
    # aspect='auto' 会自动调整单元格的宽高比以适应图形
    # cmap='viridis' 设置颜色映射方案，你可以换成 'hot', 'coolwarm', 'jet' 等
    im = plt.imshow(two_d_matrix, cmap='viridis', aspect='auto')

    # 添加颜色条 (非常重要，用于解释颜色代表的值)
    plt.colorbar(im, label='Value')

    # 添加标题和坐标轴标签
    plt.title("Heatmap from 1D Array")
    plt.xlabel("Column Index")
    plt.ylabel("Row Index")

    # 显示图形
    plt.show()