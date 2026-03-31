#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

// 如果你不需要使用 matrix 类，请将 IGNORE_MATRIX 改为 0
// #define IGNORE_MATRIX 0
#define IGNORE_MATRIX 1

#if IGNORE_MATRIX

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //       你可以任意修改 matrix 类框架的任何类成员和函数。
    //****************************

public:

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //       你可以任意修改 matrix 类框架的任何类成员和函数。
    //****************************

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // TODO: 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) : m(m_), n(n_) {
        if (m_ <= 0 || n_ <= 0) {
            throw matrix_error();
        }
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = fraction(0);
            }
        }
    }

    // TODO: 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) : m(obj.m), n(obj.n) {
        if (obj.data == nullptr) {
            data = nullptr;
            return;
        }
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }
    }

    // TODO: 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept : m(obj.m), n(obj.n), data(obj.data) {
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    // TODO: 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    // TODO: 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;

        // 释放原有内存
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }

        m = obj.m;
        n = obj.n;

        if (obj.data == nullptr) {
            data = nullptr;
            return *this;
        }

        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }

        return *this;
    }

    // TODO: 重载括号，返回矩阵的第i行(1-based)、第j列(0-based)的元素的引用。如果 i、j 不合法，抛出 matrix_error 错误。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n) {
            throw matrix_error();
        }
        return data[i-1][j];
    }

    // TODO: 重载乘号，返回矩阵乘法 lhs * rhs 的结果。如果 lhs 的列数与 rhs 的行数不相等，抛出 matrix_error 错误。
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix result(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; i++) {
            for (int j = 0; j < rhs.n; j++) {
                result.data[i][j] = fraction(0);
                for (int k = 0; k < lhs.n; k++) {
                    result.data[i][j] = result.data[i][j] + lhs.data[i][k] * rhs.data[k][j];
                }
            }
        }
        return result;
    }

    // TODO: 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
    matrix transposition() {
        if (data == nullptr || m == 0 || n == 0) {
            throw matrix_error();
        }
        matrix result(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // TODO: 返回矩阵的行列式。建议用高斯消元实现。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (data == nullptr || m == 0 || n == 0 || m != n) {
            throw matrix_error();
        }

        // 创建临时矩阵进行高斯消元
        matrix temp(*this);
        fraction det(1);

        for (int i = 0; i < m; i++) {
            // 找到主元
            int pivot = i;
            for (int j = i + 1; j < m; j++) {
                // 简单比较，找非零元素
                if (temp.data[j][i].operator==(fraction(0)) == false) {
                    if (temp.data[pivot][i].operator==(fraction(0))) {
                        pivot = j;
                    }
                }
            }

            // 如果主元为0，行列式为0
            if (temp.data[pivot][i].operator==(fraction(0))) {
                return fraction(0);
            }

            // 交换行
            if (pivot != i) {
                for (int j = 0; j < n; j++) {
                    fraction tmp = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot][j];
                    temp.data[pivot][j] = tmp;
                }
                det = det * fraction(-1);
            }

            // 消元
            for (int j = i + 1; j < m; j++) {
                fraction factor = temp.data[j][i] / temp.data[i][i];
                for (int k = i; k < n; k++) {
                    temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                }
            }

            det = det * temp.data[i][i];
        }

        return det;
    }

    // 辅助函数：获取行数和列数
    int rows() const { return m; }
    int cols() const { return n; }

    // 辅助函数：获取子矩阵（删除指定行和列）
    matrix submatrix(int row, int col) {
        if (m <= 1 || n <= 1) {
            throw matrix_error();
        }
        matrix result(m - 1, n - 1);
        int ri = 0;
        for (int i = 0; i < m; i++) {
            if (i == row) continue;
            int rj = 0;
            for (int j = 0; j < n; j++) {
                if (j == col) continue;
                result.data[ri][rj] = data[i][j];
                rj++;
            }
            ri++;
        }
        return result;
    }
};

#endif

class resistive_network {
private:


    /* 以下是建议的类成员框架，你可以选择使用或者自己实现。

    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 邻接矩阵A，电导矩阵C，Laplace矩阵(A^tCA) (具体定义见 reading_materials.pdf)
    matrix adjacency, conduction, laplace;

    */

    int interface_size, connection_size;
    matrix adjacency, conduction, laplace;

    //****************************
    // TODO: 你可以在此添加任何需要的类成员和函数。
    //****************************

public:

    //****************************
	// 注意，请勿私自修改以下函数接口的声明！
    //****************************

    // TODO: 设置电阻网络。节点数量为interface_size_，接线数量为connection_size_。
    //       对于 1<=i<=connection_size_，从节点from[i-1]到节点to[i-1]有接线，对应电阻为resistance[i-1]。
    //       保证接线使得电阻网络联通，from[i-1] < to[i-1]，resitance[i-1] > 0，均合法。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
        : interface_size(interface_size_), connection_size(connection_size_),
          adjacency(connection_size_, interface_size_),
          conduction(connection_size_, connection_size_),
          laplace(interface_size_, interface_size_) {

        // 构建邻接矩阵 A (m x n)，其中 m 是边数，n 是节点数
        // A[i][j] = 1 if edge i connects to node j+1 (outgoing), -1 if incoming, 0 otherwise
        for (int i = 0; i < connection_size_; i++) {
            for (int j = 0; j < interface_size_; j++) {
                int node_id = j + 1; // 1-based
                if (node_id == from[i]) {
                    adjacency(i+1, j) = fraction(1);
                } else if (node_id == to[i]) {
                    adjacency(i+1, j) = fraction(-1);
                } else {
                    adjacency(i+1, j) = fraction(0);
                }
            }
        }

        // 构建电导矩阵 C (m x m)，对角矩阵，C[i][i] = 1/R[i]
        for (int i = 0; i < connection_size_; i++) {
            for (int j = 0; j < connection_size_; j++) {
                if (i == j) {
                    conduction(i+1, j) = fraction(1) / resistance[i];
                } else {
                    conduction(i+1, j) = fraction(0);
                }
            }
        }

        // 构建拉普拉斯矩阵 L = A^T * C * A
        matrix At = adjacency.transposition();
        matrix CA = conduction * adjacency;
        laplace = At * CA;
    }

    ~resistive_network() = default;

    // TODO: 返回节点 interface_id1 和 interface_id2 (1-based)之间的等效电阻。
    //       保证 interface_id1 <= interface_id2 均合法。
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        // 使用公式：R_{ij} = (L_{ii}^* + L_{jj}^* - L_{ij}^* - L_{ji}^*) / det(L_n)
        // 其中 L^* 是 L 去掉最后一行和最后一列后的子矩阵的伴随矩阵元素

        if (interface_id1 == interface_id2) {
            return fraction(0);
        }

        // 去掉最后一行和最后一列，得到 L_n-1
        matrix L_reduced(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                L_reduced(i+1, j) = laplace(i+1, j);
            }
        }

        fraction det_L = L_reduced.determination();

        // 计算伴随矩阵的元素
        // L_ii^* = det(去掉第i行第i列的子矩阵) * (-1)^(i+i)
        int i1 = interface_id1 - 1; // 转为0-based
        int i2 = interface_id2 - 1;

        // 计算 L_{i1,i1}^*
        matrix sub_i1i1 = L_reduced.submatrix(i1, i1);
        fraction L_i1i1_star = sub_i1i1.determination();
        if ((i1 + i1) % 2 == 1) {
            L_i1i1_star = L_i1i1_star * fraction(-1);
        }

        // 计算 L_{i2,i2}^*
        matrix sub_i2i2 = L_reduced.submatrix(i2, i2);
        fraction L_i2i2_star = sub_i2i2.determination();
        if ((i2 + i2) % 2 == 1) {
            L_i2i2_star = L_i2i2_star * fraction(-1);
        }

        // 计算 L_{i1,i2}^* (注意：余子式的符号)
        matrix sub_i1i2 = L_reduced.submatrix(i1, i2);
        fraction L_i1i2_star = sub_i1i2.determination();
        if ((i1 + i2) % 2 == 1) {
            L_i1i2_star = L_i1i2_star * fraction(-1);
        }

        // 计算 L_{i2,i1}^*
        matrix sub_i2i1 = L_reduced.submatrix(i2, i1);
        fraction L_i2i1_star = sub_i2i1.determination();
        if ((i2 + i1) % 2 == 1) {
            L_i2i1_star = L_i2i1_star * fraction(-1);
        }

        // R_{i1,i2} = (L_{i1,i1}^* + L_{i2,i2}^* - L_{i1,i2}^* - L_{i2,i1}^*) / det(L)
        fraction numerator = L_i1i1_star + L_i2i2_star - L_i1i2_star - L_i2i1_star;
        fraction result = numerator / det_L;

        return result;
    }

    // TODO: 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size(1-based)的电压为0。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电流为 current[i-1]。
    //       保证 current 使得电阻网络有解，id < interface_size 合法。
    fraction get_voltage(int id, fraction current[]) {
        // 使用公式：U = L^(-1) * I
        // 由于 U_n = 0，我们只需要求解 L_{n-1} * U_{1..n-1} = I_{1..n-1}
        // 使用克拉默法则：U_i = det(L_i) / det(L)
        // 其中 L_i 是将 L 的第 i 列替换为 I 的矩阵

        // 去掉最后一行和最后一列
        matrix L_reduced(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                L_reduced(i+1, j) = laplace(i+1, j);
            }
        }

        fraction det_L = L_reduced.determination();

        // 构建 L_id，将第 id 列替换为电流向量
        matrix L_id(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; i++) {
            for (int j = 0; j < interface_size - 1; j++) {
                if (j == id - 1) {
                    L_id(i+1, j) = current[i];
                } else {
                    L_id(i+1, j) = L_reduced(i+1, j);
                }
            }
        }

        fraction det_L_id = L_id.determination();
        fraction voltage = det_L_id / det_L;

        return voltage;
    }


    // TODO: 在给定节点电压U的前提下，返回电阻网络的功率。
    //       对于 1<=i<=interface_size，节点i (1-based) 对应电压为 voltage[i-1]。
    //       保证 voltage 合法。
    fraction get_power(fraction voltage[]) {
        // 功率 P = U^T * L * U = sum over all edges of (U_i - U_j)^2 / R_ij
        // 更简单的方法：P = sum over all edges of (voltage difference)^2 * conductance

        fraction power(0);

        // 遍历所有边
        for (int i = 0; i < connection_size; i++) {
            // 找到这条边连接的两个节点
            int node1 = -1, node2 = -1;
            for (int j = 0; j < interface_size; j++) {
                if (adjacency(i+1, j).operator==(fraction(1))) {
                    node1 = j; // 0-based
                } else if (adjacency(i+1, j).operator==(fraction(-1))) {
                    node2 = j; // 0-based
                }
            }

            // 计算电压差
            fraction voltage_diff = voltage[node1] - voltage[node2];

            // 获取电导
            fraction conductance = conduction(i+1, i);

            // P += (U1 - U2)^2 * G
            power = power + voltage_diff * voltage_diff * conductance;
        }

        return power;
    }
};


#endif //SRC_HPP