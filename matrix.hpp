#include <array>
#include <vector>
#include <numeric>

#include <iomanip>
#include <ostream>

using namespace std;

class Matrix
{
    public:
        // Constructors
        Matrix(int M, int N): M(M), N(N), m_data(M, vector<float>(N, 0.0)) {};
        ~Matrix() = default;
        Matrix(const Matrix &) = default;
        Matrix(Matrix &&) = default;
        Matrix(const initializer_list<initializer_list<float>> &);

        // shape
        int *shape() {
            static int shape[2];
            shape[0] = M;
            shape[1] = N;
            return shape;
        };

        // assignment
        Matrix &operator=(const Matrix &) = default;
        Matrix &operator=(const initializer_list<initializer_list<float>> &);

        // index
        vector<float> &operator[](int r) {return this->row(r);};
        vector<float> &row(int r);
        vector<float> col(int c);

        // transpose
        Matrix transpose();

        // concatenate
        Matrix cat(Matrix, int);

        // print
        friend ostream &operator<<(ostream &out, Matrix);

        // addition and subtraction
        Matrix operator+(Matrix);
        void operator+=(Matrix);
        Matrix operator-(Matrix);
        void operator-=(Matrix);

        // element-wise multiplication and division
        Matrix operator*(Matrix);
        void operator*=(Matrix);
        Matrix operator/(Matrix);
        void operator/=(Matrix);

        // operation with scalar
        template <typename T>
        friend Matrix operator+(T, Matrix);
        template <typename T>
        friend Matrix operator+(Matrix, T);
        template <typename T>
        void operator+=(T);
        template <typename T>
        friend Matrix operator-(T, Matrix);
        template <typename T>
        friend Matrix operator-(Matrix, T);
        template <typename T>
        void operator-=(T);
        template <typename T>
        friend Matrix operator*(T, Matrix);
        template <typename T>
        friend Matrix operator*(Matrix, T);
        template <typename T>
        void operator*=(T);
        template <typename T>
        friend Matrix operator/(T, Matrix);
        template <typename T>
        friend Matrix operator/(Matrix, T);
        template <typename T>
        void operator/=(T);

        // matrix multiplication
        Matrix matmul(Matrix);

    private:
        int M;
        int N;
        // data
        vector<vector<float>> m_data;
};

// construction from list
Matrix::Matrix(const initializer_list<initializer_list<float>> &list) {
    M = list.size();
    N = 0;
    for (int i = 0; i < M; i++) {
        if ((list.begin()+i)->size() > N) {
            N = (list.begin()+i)->size();
        }
    }

    m_data = vector<vector<float>>(M, vector<float>(N, 0.0));
    this->operator=(list);
}

// assignment from list
Matrix &Matrix::operator=(const initializer_list<initializer_list<float>> &list) {
    int list_n_rows = list.size();
    for (int i = 0; i < M; i++) {
        int row_n_elements = (list.begin()+i)->size();
        for (int j = 0; j < N; j++) {
            if (i < list_n_rows && j < row_n_elements) {
                m_data[i][j] = *((list.begin()+i)->begin()+j);
            }
            else {
                m_data[i][j] = 0;
            }
        }
    }
    return *this;
}

// row index
vector<float> &Matrix::row(int r) {
    return m_data[r];
}

// col index
vector<float> Matrix::col(int c) {
    vector<float> col_data(this->shape()[0], 0.0);
    for (int i = 0; i < this->shape()[0]; i++) {
        col_data[i] = m_data[i][c];
    }
    return col_data;
}

// transpose
Matrix Matrix::transpose() {
    Matrix res(N, M);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            res[j][i] = m_data[i][j];
        }
    }
    return res;
}

// concatenate
Matrix Matrix::cat(Matrix other, int axis) {
    if (axis == 0) {
        int K = other.shape()[0];
        Matrix res(M+K, N);

        if (N != other.shape()[1]) {
            cout << "[!] The dimension 1 doesn't match!" << endl;
            return res;
        }

        for (int i = 0; i < M + K; i++) {
            for (int j = 0; j < N; j++) {
                if (i < M) {
                    res[i][j] = m_data[i][j];
                }
                else {
                    res[i][j] = other[i-M][j];
                }
            }
        }
        return res;
    }
    else if (axis == 1) {
        int K = other.shape()[1];
        Matrix res(M, N+K);

        if (M != other.shape()[0]) {
            cout << "[!] The dimension 0 doesn't match!" << endl;
            return res;
        }

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N + K; j++) {
                if (j < N) {
                    res[i][j] = m_data[i][j];
                }
                else {
                    res[i][j] = other[i][j-N];
                }
            }
        }
        return res;
    }
    else {
        cout << "[!] 'axis' must be 0 or 1" << endl;
        return *this;
    }
}

// print
ostream &operator<<(ostream &out, Matrix mat){
    out << "[";
    for (int i = 0; i < mat.shape()[0]; i++) {
        if (i > 0) out << " "; out << "[";
        for (int j = 0; j < mat.shape()[1]; j++) {
            if (j > 0) out << ", ";
            out << setprecision(4) << mat[i][j];
        }
        if (i < mat.shape()[0] - 1) out << "]," << endl;
        else out << "]]" << endl;
    }
    return out;
}

// addition
Matrix Matrix::operator+(Matrix other){
    Matrix res(*this);

    if (M != other.shape()[0] or N != other.shape()[1]) {
        cout << "[!] The shape of the right matrix doesn't match the left matrix!" << endl;
        return res;
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j <N; j++) {
            res[i][j] += other[i][j];
        }
    }
    return res;
}

void Matrix::operator+=(Matrix other) {
    Matrix res = this->operator+(other);
    this->operator=(res);
}

// subtraction
Matrix Matrix::operator-(Matrix other) {
    Matrix res(*this);

    if (M != other.shape()[0] or N != other.shape()[1]) {
        cout << "[!] The shape of the right matrix doesn't match the left matrix!" << endl;
        return res;
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            res[i][j] -= other[i][j];
        }
    }
    return res;
}

void Matrix::operator-=(Matrix other) {
    Matrix res = this->operator-(other);
    this->operator=(res);
}

// element-wise multiplication
Matrix Matrix::operator*(Matrix other) {
    Matrix res(*this);

    if (M != other.shape()[0] or N != other.shape()[1]) {
        cout << "[!] The shape of the right matrix doesn't match the left matrix!" << endl;
        return res;
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            res[i][j] *= other[i][j];
        }
    }
    return res;
}

void Matrix::operator*=(Matrix other) {
    Matrix res = this->operator*(other);
    this->operator=(res);
}

// element-wise division
Matrix Matrix::operator/(Matrix other) {
    Matrix res(*this);

    if (M != other.shape()[0] or N != other.shape()[1]) {
        cout << "[!] The shape of the right matrix doesn't match the left matrix!" << endl;
        return res;
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            res[i][j] /= other[i][j];
        }
    }
    return res;
}

void Matrix::operator/=(Matrix other) {
    Matrix res = this->operator/(other);
    this->operator=(res);
}

// operator with scalar, scalar add matrix
template <typename T>
Matrix operator+(T val, Matrix mat) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] = val + mat[i][j];
        }
    }
    return res;
}

// operator with scalar, matrix add scalar
template <typename T>
Matrix operator+(Matrix mat, T val) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] += val;
        }
    }
    return res;
}

template <typename T>
void Matrix::operator+=(T val) {
    Matrix res = *this + val;
    this->operator=(res);
}

// operator with scalar, scalar subtract matrix
template <typename T>
Matrix operator-(T val, Matrix mat) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] = val - mat[i][j];
        }
    }
    return res;
}

// operator with scalar, matrix subtract scalar
template <typename T>
Matrix operator-(Matrix mat, T val) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] -= val;
        }
    }
    return res;
}

template <typename T>
void Matrix::operator-=(T val) {
    Matrix res = *this - val;
    this->operator=(res);
}

// operator with scalar, scalar multiply matrix
template <typename T>
Matrix operator*(T val, Matrix mat) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] = val * mat[i][j];
        }
    }
    return res;
}

// operator with scalar, matrix multiply scalar
template <typename T>
Matrix operator*(Matrix mat, T val) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] *= val;
        }
    }
    return res;
}

template <typename T>
void Matrix::operator*=(T val) {
    Matrix res = *this * val;
    this->operator=(res);
}

// operator with scalar, scalar divide matrix
template <typename T>
Matrix operator/(T val, Matrix mat) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] = val / mat[i][j];
        }
    }
    return res;
}

// operator with scalar, matrix divide scalar
template <typename T>
Matrix operator/(Matrix mat, T val) {
    Matrix res(mat);
    for (int i = 0; i < res.shape()[0]; i++) {
        for (int j = 0; j < res.shape()[1]; j++) {
            res[i][j] /= val;
        }
    }
    return res;
}

template <typename T>
void Matrix::operator/=(T val) {
    Matrix res = *this / val;
    this->operator=(res);
}

// matrix multiplication, [M x N] * [N x K] = [M x K]
Matrix Matrix::matmul(Matrix other) {
    int K = other.shape()[1];
    Matrix res(M, K);

    if (N != other.shape()[0]) {
        cout << "[!] The shape of the right matrix doesn't match the left matrix!" << endl;
        return res;
    }
    
    vector<float> row;
    vector<float> col;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < K; j++) {
            row = this->row(i);
            col = other.col(j);
            res[i][j] = inner_product(row.begin(), row.end(), col.begin(), 0.0);
        }
    }
    return res;
}
