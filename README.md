#### Brief

A simple C++ implementation of matrix operations.



#### Usage Example

```c++
#include <iostream>
#include "matrix.hpp"

using namespace std;

int main(){
    // init a matrix from a 2D list
    Matrix a = {{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}};
    Matrix b = {{4.11, 6.06}, {8.07, 7.20}};

    // init a 2x2 zero matrix
    Matrix c(2, 2);

    // print matrix to the console
    cout << a << endl;
    cout << b << endl;
    cout << c << endl;

    // index
    float val = b[0][0];

    // transpose
    b = b.transpose();

    // concatenate
    Matrix x = b.cat(c, 0);     // concatenate matrices along axis 0
    Matrix y = a.cat(b, 1);     // concatenate matrices along axis 1

    Matrix sum = b + c;         // element-wise addition
    Matrix diff = b - c;        // element-wise subtraction
    Matrix prod = b * c;        // element-wise multiplication
    Matrix quot = b / b;        // element-wise division

    // compound assignment operator
    b += c;
    b -= c;
    b *= c;
    a /= a;

    // operations with scalar;
    a = a + 1;
    b = 4.11 - b;
    c *= 6;
    x = y / 2.0;

    // matrix multiplication
    Matrix mat_prod = b.matmul(c);

    return 0;
}
```
