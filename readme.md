# Numerical Methods Console Application

## 👋 Welcome to FindX

**FindX** is a Numerical Methods laboratory project for the course **CSE 2208 — Numerical Method**. It organizes the numerical techniques taught in the course into a clear and structured repository, combining **mathematical theory**, **algorithmic steps**, and **C++ implementations**. The project is designed for academic learning, and quick revision, with a focus on clarity, correctness, and practical application.

---

## Project Structure

```text
Numerical Methods
├── readme.md
│
├── A. Solution of Linear Equations
│   ├── Gauss Elimination
│   ├── Gauss Jordan
│   ├── LU Decomposition
│   ├── Matrix Inverse
│   └── Iterative Methods
│       ├── Jacobi
│       └── GaussSeidel
│
├── B. Solution of Non-Linear Equations
│   ├── Bisection
│   ├── False_Position
│   ├── Secant
│   └── Newton_Raphson
│
├── C. Interpolation and Approximation
│   ├── Newton Forward
│   └── Newton Backward
│
├── D. Numerical Differentiation/
│   ├── Newton Forward Differentiation/
│
├── E. Solution of Differential Equations/
│   └── Runge Kutta/
│
├── F. Numerical Integration
│   ├── Simpson 1 by 3
│   └── Simpson 3 by 8
│
└── G. Curve Fitting
    ├── Linear
    ├── Polynomial
    └── Transcendental
```

---

## Environment Setup

To run the numerical method implementations in this repository, you will need a standard C++ development environment. All code is written in **Standard C++ (C++11 or later)** to ensure compatibility.

### Prerequisites
* **C++ Compiler:** G++ (via MinGW for Windows, or GCC for Linux/macOS).
* **Code Editor/IDE:** Visual Studio Code (recommended), CLion, or Dev-C++.

You can compile and run the programs using any IDE or directly from the terminal.

**Using Terminal (Command Line):**
1.  Navigate to the directory of the specific method you want to run.
2.  Compile the source code:
    ```bash
    g++ main.cpp -o run
    ```
3.  Execute the program:
    * **Windows:** `run.exe`
    * **Linux/Mac:** `./run`

**Using VS Code:**
* Open the folder containing the source code.
* Ensure the C/C++ extension is installed.
* Press `Ctrl + Shift + B` (Build) or use the "Run Code" button if you have Code Runner installed.
---

## Table of Contents

### A. Solution of Linear Equations
1.  [Gauss Elimination Method](#gauss-elimination)
2.  [Gauss Jordan Method](#gauss-jordan)
3.  [LU Decomposition Method](#lu-decomposition)
4.  [Jacobi Method](#jacobi)
5.  [GaussSeidel Method](#gaussseidel)

### B. Solution of Non-Linear Equations
1.  [Bisection Method](#bisection)
2.  [False Position Method](#false_position)
3.  [Secant Method](#secant)
4.  [Newton Raphson Method](#newton_raphson)

### C. Interpolation and Approximation
1.  [Newton Forward Interpolation](#newton-forward)
2.  [Newton Backward Interpolation](#newton-backward)

### D. Numerical Differentiation
1.  [Newton Forward Differentiation](#newton-forward-differentiation)


### E. Solution of Differential Equations
1.  [Runge Kutta Method](#runge-kutta)

### F. Numerical Integration
1.  [Simpson 1 by 3](#simpson-1-by-3)
2.  [Simpson 3 by 8](#simpson-3-by-8)

### G. Curve Fitting Methods
1.  [Linear Equation](#linear)
2.  [Polynomial Equation](#polynomial)
3.  [Transcendental Equation](#transcendental)

---



## A. Solution of Linear Equations

<a id="gauss-elimination"></a>
### 1. Gauss Elimination Method

**Theory**
Gauss Elimination is a direct method that converts a system of linear equations into an upper triangular system using forward elimination. After the matrix becomes upper triangular, back substitution is applied to determine the values of the unknowns. Gauss Elimination produces one of the following outcomes:
1. Unique Solution
2. No Solution (Inconsistent System)
3. Infinitely Many Solutions (Dependent System)

The process of finding the solution includes writing the augmented matrix [A|B]. Applying forward elimination to make all elements below the diagonal equal to zero. After forward elimination, the matrix becomes upper triangular. Applying back substitution starting from the last equation to find the solutions.

Detecting:
* If a row becomes [0 0 0 | c] with c != 0 => No Solution
* If a row becomes [0 0 0 | 0] => Infinite Solutions

**Algorithm**
Input:
    Number of equations n
    Augmented matrix A[n][n+1]

Steps:
1. Forward Elimination:
   For each column k = 1 to n-1:
    a. Pivoting:
        Finding row with maximum absolute value in column k
        Swapping the pivot row with row k
    b. Elimination:
        For each row i = k+1 to n:
        m = A[i][k] / A[k][k]
        For each column j = k to n+1:
        A[i][j] = A[i][j] - m * A[k][j]

2. Checking for No Solution:
   If a row is [0 0 0 | c] where c != 0 => No solution

3. Checking for Infinite Solutions:
   If a row is [0 0 0 ... 0 | 0] and rank(A) = rank(A|b) < n => Infinite solutions

4. Back Substitution:
   For i = n down to 1:
   sum = Σ (A[i][j] * x[j]) for j = i+1 to n
   x[i] = (A[i][n+1] - sum) / A[i][i]

**Pseudocode**
```text
START
Input n
Input augmented matrix A[n][n+1]

# FORWARD ELIMINATION
for k = 1 to n-1:
    # Pivoting
    pivot_row = k
    for i = k+1 to n:
        if abs(A[i][k]) > abs(A[pivot_row][k]):
            pivot_row = i
    swap rows A[k] and A[pivot_row]

    if A[k][k] == 0:
        continue

    # Eliminate rows below pivot
    for i = k+1 to n:
        m = A[i][k] / A[k][k]
        for j = k to n+1:
            A[i][j] = A[i][j] - m * A[k][j]

# CHECK SOLUTION TYPE
infinite_flag = false
no_solution_flag = false

for i = 1 to n:
    if all A[i][1..n] == 0 AND A[i][n+1] != 0:
        no_solution_flag = true
    if all A[i][1..n] == 0 AND A[i][n+1] == 0:
        infinite_flag = true

if no_solution_flag:
    PRINT "No Solution"
    STOP

if infinite_flag:
    PRINT "Infinite Solutions"
    STOP

# BACK SUBSTITUTION
Create vector x[n]

for i = n to 1 step -1:
    sum = 0
    for j = i+1 to n:
        sum = sum + A[i][j] * x[j]
    x[i] = (A[i][n+1] - sum) / A[i][i]

PRINT "Unique Solution: ", x
STOP
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

void print_matrix(vector<vector<double>> &a, ostream &out)
{
    int n = a.size();
    int m = a[0].size();
    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < m; j++)
        {
            out << fixed << setprecision(6) << a[i][j] << " ";
        }
        out << "\n";
    }
    out << "\n";
}

int gauss_elimination(vector<vector<double>> a, vector<double> &root, ostream &out)
{
    int n = a.size();

    out << "Initial Augmented Matrix:\n";
    print_matrix(a, out);

    for (int i = 0; i < n; i++)
    {
        int pivot = i;

        for (int r = i + 1; r < n; r++)
            if (fabs(a[r][i]) > fabs(a[pivot][i]))
                pivot = r;

        if (fabs(a[pivot][i]) < 1e-12)
            continue;

        swap(a[i], a[pivot]);

        for (int r = i + 1; r < n; r++)
        {
            if (fabs(a[r][i]) < 1e-12)
                continue;

            double factor = a[r][i] / a[i][i];
            for (int c = i; c <= n; c++)
                a[r][c] -= factor * a[i][c];
        }

        out << "After eliminating column " << i + 1 << ":\n";
        print_matrix(a, out);
    }

    for (int i = 0; i < n; i++)
    {
        bool zeroRow = true;
        for (int j = 0; j < n; j++)
            if (fabs(a[i][j]) > 1e-12)
                zeroRow = false;

        if (zeroRow && fabs(a[i][n]) > 1e-12)
            return -1;
        if (zeroRow && fabs(a[i][n]) < 1e-12)
            return 0;
    }

    root.assign(n, 0.0);
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = a[i][n];
        for (int j = i + 1; j < n; j++)
            sum -= a[i][j] * root[j];

        if (fabs(a[i][i]) < 1e-12)
            return 0;
        root[i] = sum / a[i][i];
    }

    return 1;
}

int main()
{
    ifstream fin("inputgauss.txt");
    ofstream fout("outputgauss.txt");

    int T;
    fin >> T;

    for (int t = 1; t <= T; t++)
    {
        int n;
        fin >> n;

        vector<vector<double>> a(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> a[i][j];

        fout << t << " (Gauss Elimination)\n";

        vector<double> root;
        int status = gauss_elimination(a, root, fout);

        if (status == 1)
        {
            fout << "Unique Solution:\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << root[i] << "\n";
        }
        else if (status == 0)
            fout << "Infinite Solutions\n";
        else
            fout << "No Solution\n";

        fout << "\n\n";
    }

    return 0;
}

```

#### Input
```
9
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
3
1 1 1 6
0 2 5 -4
2 5 -1 27
3
4 -2 1 9
1 1 1 6
3 4 -2 -1
3
1 -2 3 4
2 -4 6 10
3 1 -1 2
3
1 1 1 3
2 2 2 7
1 -1 2 4
3
3 -6 9 12
1 -2 3 5
2 4 -1 7
3
1 2 -1 3
2 4 -2 6
3 6 -3 9
3
2 1 -3 4
4 2 -6 8
1 0.5 -1.5 2
3
1 -1 2 3
2 -2 4 6
3 -3 6 9

```

#### Output
```
System 1 (Gauss-Jordan)
Initial Augmented Matrix:
2.000000 1.000000 -1.000000 8.000000 
-3.000000 -1.000000 2.000000 -11.000000 
-2.000000 1.000000 2.000000 -3.000000 

After pivot normalization in row 1:
1.000000 0.333333 -0.666667 3.666667 
0.000000 0.333333 0.333333 0.666667 
0.000000 1.666667 0.666667 4.333333 

After pivot normalization in row 2:
1.000000 0.000000 -0.800000 2.800000 
0.000000 1.000000 0.400000 2.600000 
0.000000 0.000000 0.200000 -0.200000 

After pivot normalization in row 3:
1.000000 0.000000 0.000000 2.000000 
0.000000 1.000000 0.000000 3.000000 
0.000000 0.000000 1.000000 -1.000000 

Unique Solution:
x1 = 2.000000
x2 = 3.000000
x3 = -1.000000


System 2 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 1.000000 1.000000 6.000000 
0.000000 2.000000 5.000000 -4.000000 
2.000000 5.000000 -1.000000 27.000000 

After pivot normalization in row 1:
1.000000 2.500000 -0.500000 13.500000 
0.000000 2.000000 5.000000 -4.000000 
0.000000 -1.500000 1.500000 -7.500000 

After pivot normalization in row 2:
1.000000 0.000000 -6.750000 18.500000 
0.000000 1.000000 2.500000 -2.000000 
0.000000 0.000000 5.250000 -10.500000 

After pivot normalization in row 3:
1.000000 0.000000 0.000000 5.000000 
0.000000 1.000000 0.000000 3.000000 
0.000000 0.000000 1.000000 -2.000000 

Unique Solution:
x1 = 5.000000
x2 = 3.000000
x3 = -2.000000


System 3 (Gauss-Jordan)
Initial Augmented Matrix:
4.000000 -2.000000 1.000000 9.000000 
1.000000 1.000000 1.000000 6.000000 
3.000000 4.000000 -2.000000 -1.000000 

After pivot normalization in row 1:
1.000000 -0.500000 0.250000 2.250000 
0.000000 1.500000 0.750000 3.750000 
0.000000 5.500000 -2.750000 -7.750000 

After pivot normalization in row 2:
1.000000 0.000000 0.000000 1.545455 
0.000000 1.000000 -0.500000 -1.409091 
0.000000 0.000000 1.500000 5.863636 

After pivot normalization in row 3:
1.000000 0.000000 0.000000 1.545455 
0.000000 1.000000 0.000000 0.545455 
0.000000 0.000000 1.000000 3.909091 

Unique Solution:
x1 = 1.545455
x2 = 0.545455
x3 = 3.909091


System 4 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 -2.000000 3.000000 4.000000 
2.000000 -4.000000 6.000000 10.000000 
3.000000 1.000000 -1.000000 2.000000 

After pivot normalization in row 1:
1.000000 0.333333 -0.333333 0.666667 
0.000000 -4.666667 6.666667 8.666667 
0.000000 -2.333333 3.333333 3.333333 

After pivot normalization in row 2:
1.000000 0.000000 0.142857 1.285714 
-0.000000 1.000000 -1.428571 -1.857143 
0.000000 0.000000 0.000000 -1.000000 

No Solution


System 5 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 1.000000 1.000000 3.000000 
2.000000 2.000000 2.000000 7.000000 
1.000000 -1.000000 2.000000 4.000000 

After pivot normalization in row 1:
1.000000 1.000000 1.000000 3.500000 
0.000000 0.000000 0.000000 -0.500000 
0.000000 -2.000000 1.000000 0.500000 

After pivot normalization in row 2:
1.000000 0.000000 1.500000 3.750000 
-0.000000 1.000000 -0.500000 -0.250000 
0.000000 0.000000 0.000000 -0.500000 

No Solution


System 6 (Gauss-Jordan)
Initial Augmented Matrix:
3.000000 -6.000000 9.000000 12.000000 
1.000000 -2.000000 3.000000 5.000000 
2.000000 4.000000 -1.000000 7.000000 

After pivot normalization in row 1:
1.000000 -2.000000 3.000000 4.000000 
0.000000 0.000000 0.000000 1.000000 
0.000000 8.000000 -7.000000 -1.000000 

After pivot normalization in row 2:
1.000000 0.000000 1.250000 3.750000 
0.000000 1.000000 -0.875000 -0.125000 
0.000000 0.000000 0.000000 1.000000 

No Solution


System 7 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 2.000000 -1.000000 3.000000 
2.000000 4.000000 -2.000000 6.000000 
3.000000 6.000000 -3.000000 9.000000 

After pivot normalization in row 1:
1.000000 2.000000 -1.000000 3.000000 
0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 

Infinite Solutions


System 8 (Gauss-Jordan)
Initial Augmented Matrix:
2.000000 1.000000 -3.000000 4.000000 
4.000000 2.000000 -6.000000 8.000000 
1.000000 0.500000 -1.500000 2.000000 

After pivot normalization in row 1:
1.000000 0.500000 -1.500000 2.000000 
0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 

Infinite Solutions


System 9 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 -1.000000 2.000000 3.000000 
2.000000 -2.000000 4.000000 6.000000 
3.000000 -3.000000 6.000000 9.000000 

After pivot normalization in row 1:
1.000000 -1.000000 2.000000 3.000000 
0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 

Infinite Solutions



```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/A.%20Solution%20of%20Linear%20Equations/Gauss%20Elimination)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Gauss Elimination - Wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination)
- [Gauss Elimination Method - GeeksforGeeks](https://www.geeksforgeeks.org/gaussian-elimination/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="gauss-jordan"></a>
### 2. Gauss Jordan Method

**Theory**
Gauss-Jordan Elimination is an extended form of Gauss Elimination. Instead of producing an upper triangular matrix, it reduces the augmented matrix directly to reduced row echelon form (RREF). This eliminates the need for back substitution.

The method for solving system of linear equations include writing the augmented matrix [A|B]. Converting the matrix to upper triangular form (like Gauss elimination). Continuing eliminating values above the diagonal to create a diagonal matrix. Converting diagonal elements to 1 by dividing the row.
The matrix becomes:
    [ I | X ]
Where I = identity matrix, X = solutions.

**Algorithm**
Input:
    Number of equations n
    Augmented matrix A[n][n+1]

Steps:
1. For each pivot column k = 1 to n:
    a. Pivoting:
        Find row with maximum absolute value in column k
        Swap it with row k
    b. Normalize Pivot:
        pivot = A[k][k]
        Divide entire row k by pivot to make pivot = 1
    c. Eliminate All Other Rows:
        For each row i != k:
            m = A[i][k]
            For each column j = k to n+1:
            A[i][j] = A[i][j] - m * A[k][j]

2. Detect Solution Type:
    No solution if row is [0 0 0 | c], c != 0
    Infinite solutions if row is [0 0 0 | 0] and rank < n

3. Extract solution:
    If unique, solution x[i] = A[i][n+1]

**Pseudocode**
```text
Algorithm GaussJordan(AugmentedMatrix M):
Input : M is an m x (n+1) augmented matrix [A | b]
Output: RREF of M and type of solution (unique / none / infinite)

# Step 1 - Initialization
pivot_row <- 0
pivot_col <- 0

# Step 2 - Looping over all columns except the last one (which is b)
while pivot_row < m AND pivot_col < n:

    # Step 2.1 - Finding a row with a non-zero entry in pivot_col
    row_with_pivot <- -1
    for r from pivot_row to m-1:
        if M[r][pivot_col] != 0:
            row_with_pivot <- r
            break
    
    # If no pivot found in this column => move to next column
    if row_with_pivot == -1:
        pivot_col <- pivot_col + 1
        continue

    # Step 2.2 - Swap pivot row into correct position
    SwapRows(M, pivot_row, row_with_pivot)

    # Step 2.3 - Normalize pivot row (make pivot = 1)
    pivot_value <- M[pivot_row][pivot_col]
    for c from pivot_col to n:
        M[pivot_row][c] <- M[pivot_row][c] / pivot_value

    # Step 2.4 - Eliminate all other rows
    for r from 0 to m-1:
        if r != pivot_row:
            factor <- M[r][pivot_col]
            for c from pivot_col to n:
                M[r][c] <- M[r][c] - factor * M[pivot_row][c]

    # Move to next pivot location
    pivot_row <- pivot_row + 1
    pivot_col <- pivot_col + 1

# Step 3 - Analyze RREF for solutions
contradiction <- false
for r from 0 to m-1:
    if (all M[r][0..n-1] == 0) AND (M[r][n] != 0):
        contradiction <- true

if contradiction:
    return ("No solution", M)

# Count pivot rows
rank <- number of rows that contain a pivot (leading 1)

if rank = n:
    return ("Unique solution", M)
else:
    return ("Infinite solutions", M)
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

void print_matrix(vector<vector<double>> &a, ostream &out)
{
    int n = a.size();
    int m = a[0].size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            out << fixed << setprecision(6) << a[i][j] << " ";
        }
        out << "\n";
    }
    out << "\n";
}

int gauss_jordan(vector<vector<double>> a, vector<double> &root, ostream &out)
{
    int n = a.size();

    out << "Initial Augmented Matrix:\n";
    print_matrix(a, out);

    for (int i = 0; i < n; i++)
    {
        int pivot = i;

        for (int r = i + 1; r < n; r++)
            if (fabs(a[r][i]) > fabs(a[pivot][i]))
                pivot = r;

        if (fabs(a[pivot][i]) < 1e-12)
            continue;

        swap(a[i], a[pivot]);

        double div = a[i][i];
        for (int c = 0; c <= n; c++)
            a[i][c] /= div;

        for (int r = 0; r < n; r++)
        {
            if (r == i)
                continue;
            double factor = a[r][i];
            for (int c = 0; c <= n; c++)
                a[r][c] -= factor * a[i][c];
        }

        out << "After pivot normalization in row " << i + 1 << ":\n";
        print_matrix(a, out);
    }

    for (int i = 0; i < n; i++)
    {
        bool zeroRow = true;
        for (int j = 0; j < n; j++)
            if (fabs(a[i][j]) > 1e-12)
                zeroRow = false;

        if (zeroRow && fabs(a[i][n]) > 1e-12)
            return -1;
        if (zeroRow && fabs(a[i][n]) < 1e-12)
            return 0;
    }

    root.assign(n, 0.0);
    for (int i = 0; i < n; i++)
        root[i] = a[i][n];

    return 1;
}

int main()
{
    ifstream fin("inputgauss.txt");
    ofstream fout("outputgauss.txt");

    int T;
    fin >> T;

    for (int t = 1; t <= T; t++)
    {
        int n;
        fin >> n;

        vector<vector<double>> a(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> a[i][j];

        fout << "System " << t << " (Gauss-Jordan)\n";

        vector<double> root;
        int status = gauss_jordan(a, root, fout);

        if (status == 1)
        {
            fout << "Unique Solution:\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << root[i] << "\n";
        }
        else if (status == 0)
            fout << "Infinite Solutions\n";
        else
            fout << "No Solution\n";

        fout << "\n\n";
    }

    return 0;
}

```

#### Input
```
9
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
3
1 1 1 6
0 2 5 -4
2 5 -1 27
3
4 -2 1 9
1 1 1 6
3 4 -2 -1
3
1 -2 3 4
2 -4 6 10
3 1 -1 2
3
1 1 1 3
2 2 2 7
1 -1 2 4
3
3 -6 9 12
1 -2 3 5
2 4 -1 7
3
1 2 -1 3
2 4 -2 6
3 6 -3 9
3
2 1 -3 4
4 2 -6 8
1 0.5 -1.5 2
3
1 -1 2 3
2 -2 4 6
3 -3 6 9

```

#### Output
```
System 1 (Gauss-Jordan)
Initial Augmented Matrix:
2.000000 1.000000 -1.000000 8.000000 
-3.000000 -1.000000 2.000000 -11.000000 
-2.000000 1.000000 2.000000 -3.000000 

After pivot normalization in row 1:
1.000000 0.333333 -0.666667 3.666667 
0.000000 0.333333 0.333333 0.666667 
0.000000 1.666667 0.666667 4.333333 

After pivot normalization in row 2:
1.000000 0.000000 -0.800000 2.800000 
0.000000 1.000000 0.400000 2.600000 
0.000000 0.000000 0.200000 -0.200000 

After pivot normalization in row 3:
1.000000 0.000000 0.000000 2.000000 
0.000000 1.000000 0.000000 3.000000 
0.000000 0.000000 1.000000 -1.000000 

Unique Solution:
x1 = 2.000000
x2 = 3.000000
x3 = -1.000000


System 2 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 1.000000 1.000000 6.000000 
0.000000 2.000000 5.000000 -4.000000 
2.000000 5.000000 -1.000000 27.000000 

After pivot normalization in row 1:
1.000000 2.500000 -0.500000 13.500000 
0.000000 2.000000 5.000000 -4.000000 
0.000000 -1.500000 1.500000 -7.500000 

After pivot normalization in row 2:
1.000000 0.000000 -6.750000 18.500000 
0.000000 1.000000 2.500000 -2.000000 
0.000000 0.000000 5.250000 -10.500000 

After pivot normalization in row 3:
1.000000 0.000000 0.000000 5.000000 
0.000000 1.000000 0.000000 3.000000 
0.000000 0.000000 1.000000 -2.000000 

Unique Solution:
x1 = 5.000000
x2 = 3.000000
x3 = -2.000000


System 3 (Gauss-Jordan)
Initial Augmented Matrix:
4.000000 -2.000000 1.000000 9.000000 
1.000000 1.000000 1.000000 6.000000 
3.000000 4.000000 -2.000000 -1.000000 

After pivot normalization in row 1:
1.000000 -0.500000 0.250000 2.250000 
0.000000 1.500000 0.750000 3.750000 
0.000000 5.500000 -2.750000 -7.750000 

After pivot normalization in row 2:
1.000000 0.000000 0.000000 1.545455 
0.000000 1.000000 -0.500000 -1.409091 
0.000000 0.000000 1.500000 5.863636 

After pivot normalization in row 3:
1.000000 0.000000 0.000000 1.545455 
0.000000 1.000000 0.000000 0.545455 
0.000000 0.000000 1.000000 3.909091 

Unique Solution:
x1 = 1.545455
x2 = 0.545455
x3 = 3.909091


System 4 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 -2.000000 3.000000 4.000000 
2.000000 -4.000000 6.000000 10.000000 
3.000000 1.000000 -1.000000 2.000000 

After pivot normalization in row 1:
1.000000 0.333333 -0.333333 0.666667 
0.000000 -4.666667 6.666667 8.666667 
0.000000 -2.333333 3.333333 3.333333 

After pivot normalization in row 2:
1.000000 0.000000 0.142857 1.285714 
-0.000000 1.000000 -1.428571 -1.857143 
0.000000 0.000000 0.000000 -1.000000 

No Solution


System 5 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 1.000000 1.000000 3.000000 
2.000000 2.000000 2.000000 7.000000 
1.000000 -1.000000 2.000000 4.000000 

After pivot normalization in row 1:
1.000000 1.000000 1.000000 3.500000 
0.000000 0.000000 0.000000 -0.500000 
0.000000 -2.000000 1.000000 0.500000 

After pivot normalization in row 2:
1.000000 0.000000 1.500000 3.750000 
-0.000000 1.000000 -0.500000 -0.250000 
0.000000 0.000000 0.000000 -0.500000 

No Solution


System 6 (Gauss-Jordan)
Initial Augmented Matrix:
3.000000 -6.000000 9.000000 12.000000 
1.000000 -2.000000 3.000000 5.000000 
2.000000 4.000000 -1.000000 7.000000 

After pivot normalization in row 1:
1.000000 -2.000000 3.000000 4.000000 
0.000000 0.000000 0.000000 1.000000 
0.000000 8.000000 -7.000000 -1.000000 

After pivot normalization in row 2:
1.000000 0.000000 1.250000 3.750000 
0.000000 1.000000 -0.875000 -0.125000 
0.000000 0.000000 0.000000 1.000000 

No Solution


System 7 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 2.000000 -1.000000 3.000000 
2.000000 4.000000 -2.000000 6.000000 
3.000000 6.000000 -3.000000 9.000000 

After pivot normalization in row 1:
1.000000 2.000000 -1.000000 3.000000 
0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 

Infinite Solutions


System 8 (Gauss-Jordan)
Initial Augmented Matrix:
2.000000 1.000000 -3.000000 4.000000 
4.000000 2.000000 -6.000000 8.000000 
1.000000 0.500000 -1.500000 2.000000 

After pivot normalization in row 1:
1.000000 0.500000 -1.500000 2.000000 
0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 

Infinite Solutions


System 9 (Gauss-Jordan)
Initial Augmented Matrix:
1.000000 -1.000000 2.000000 3.000000 
2.000000 -2.000000 4.000000 6.000000 
3.000000 -3.000000 6.000000 9.000000 

After pivot normalization in row 1:
1.000000 -1.000000 2.000000 3.000000 
0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 

Infinite Solutions



```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/A.%20Solution%20of%20Linear%20Equations/Gauss%20Jordan)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Gauss-Jordan Elimination - Wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination#Gauss%E2%80%93Jordan_elimination)
- [Gauss Jordan Method - GeeksforGeeks](https://www.geeksforgeeks.org/gauss-jordan-elimination-method/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="lu-decomposition"></a>
### 3. LU Decomposition Method 

**Theory**
1. **The main idea**: Given a square matrix (**A**), it will be rewritten as $A = LU$, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
2. **Using LU to solve AX = B**:
    Once we have $A = LU$, we can write $LUx = b$.
    We introduce a helper vector **Y**: $Ly = b$ and $Ux = y$.
    First, solve $Ly = b$ using **forward substitution**.
    Then, solve $Ux = y$ using **backward substitution**.
3. **Why pivoting is often needed**:
   Sometimes LU decomposition runs into trouble if a diagonal element (called a **pivot**) becomes zero. To fix this, we often swap rows to bring a better pivot into position. $PA = LU$. This is called **partial pivoting**.

**Algorithm**
For each column/step ($k = 1$) to ($n$):
1.  Compute the ($k$)-th row of (**U**)
2.  Compute the ($k$)-th column of (**L**)

Mathematically:
$$U_{k,j} = A_{k,j} - \sum_{s=1}^{k-1} L_{k,s} U_{s,j} \quad , \quad j = k, \dots, n$$
$$L_{i,k} = \frac{A_{i,k} - \sum_{s=1}^{k-1} L_{i,s} U_{s,k}}{U_{k,k}} \quad , \quad i = k+1, \dots, n$$

**Pseudocode**
```cpp
function solve_LU(Matrix A, Vector b):
    n = size(A)
    L = IdentityMatrix(n)
    U = ZeroMatrix(n)
    
    for i from 0 to n-1:
        pivot_row = find_max_in_column(i, from: i to n-1)
        swap_rows(A, i, pivot_row)
        swap_rows(b, i, pivot_row) 
        
        for k from i to n-1:
            prod_sum = 0
            for j from 0 to i-1:
                prod_sum += L[i][j] * U[j][k]
            U[i][k] = A[i][k] - prod_sum
            
        for k from i to n-1:
            if i == k:
                L[i][i] = 1
            else:
                prod_sum = 0
                for j from 0 to i-1:
                    prod_sum += L[k][j] * U[j][i]
                L[k][i] = (A[k][i] - prod_sum) / U[i][i]

    for i from 0 to n-1:
        if abs(U[i][i]) < epsilon:
            return "No unique solution"

    Vector y(n)
    for i from 0 to n-1:
        sum = 0
        for j from 0 to i-1:
            sum += L[i][j] * y[j]
        y[i] = b[i] - sum

    Vector x(n)
    for i from n-1 down to 0:
        sum = 0
        for j from i+1 to n-1:
            sum += U[i][j] * x[j]
        x[i] = (y[i] - sum) / U[i][i]
        
    return x
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream infile("inputlu.txt");
    ofstream outfile("outputlu.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                infile >> a[i][j];
            infile >> b[i];
        }

        vector<vector<double>> l(n, vector<double>(n, 0));
        vector<vector<double>> u(n, vector<double>(n, 0));

        for (int i = 0; i < n; i++)
        {
            for (int k = i; k < n; k++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += l[i][j] * u[j][k];
                u[i][k] = a[i][k] - sum;
            }

            for (int k = i; k < n; k++)
            {
                if (i == k)
                    l[i][i] = 1;
                else
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += l[k][j] * u[j][i];

                    if (abs(u[i][i]) > 1e-9)
                        l[k][i] = (a[k][i] - sum) / u[i][i];
                    else
                        l[k][i] = 0;
                }
            }
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Lower Triangular Matrix :" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                outfile << fixed << setprecision(3) << setw(8) << l[i][j] << " ";
            outfile << endl;
        }
        outfile << endl;

        outfile << "Upper Triangular Matrix :" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                outfile << fixed << setprecision(3) << setw(8) << u[i][j] << " ";
            outfile << endl;
        }

        vector<double> y(n);
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += l[i][j] * y[j];
            y[i] = b[i] - sum;
        }

        bool unique = true;
        bool inconsistent = false;

        for (int i = 0; i < n; i++)
        {
            if (abs(u[i][i]) < 1e-9)
            {
                unique = false;
                bool allZero = true;
                for (int j = i + 1; j < n; j++)
                {
                    if (abs(u[i][j]) > 1e-9)
                    {
                        allZero = false;
                        break;
                    }
                }

                if (allZero && abs(y[i]) > 1e-9)
                {
                    inconsistent = true;
                    break;
                }
            }
        }

        if (unique)
        {
            vector<double> x(n);
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i + 1; j < n; j++)
                    sum += u[i][j] * x[j];
                x[i] = (y[i] - sum) / u[i][i];
            }

            outfile << "\nSolution:" << endl;
            for (int i = 0; i < n; i++)
                outfile << "x" << i + 1 << " = " << fixed << setprecision(4) << x[i] << endl;

            outfile << "\nThe system has a unique solution." << endl;
        }
        else
        {
            if (inconsistent)
                outfile << "\nThe system has no solution." << endl;
            else
                outfile << "\nThe system has infinite solutions." << endl;
        }
        outfile << endl
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
3
3 -0.1 -0.2 7.85
0.1 7 -0.3 -19.3
0.3 -0.2 10 71.4
3
1 1 1 6
1 1 1 7
1 2 3 14
3
1 1 1 6
2 2 2 12
1 2 3 14

```

#### Output
```
Case 1:
Lower Triangular Matrix :
   1.000    0.000    0.000 
   0.033    1.000    0.000 
   0.100   -0.027    1.000 

Upper Triangular Matrix :
   3.000   -0.100   -0.200 
   0.000    7.003   -0.293 
   0.000    0.000   10.012 

Solution:
x1 = 3.0000
x2 = -2.5000
x3 = 7.0000

The system has a unique solution.


Case 2:
Lower Triangular Matrix :
   1.000    0.000    0.000 
   1.000    1.000    0.000 
   1.000    0.000    1.000 

Upper Triangular Matrix :
   1.000    1.000    1.000 
   0.000    0.000    0.000 
   0.000    0.000    2.000 

The system has no solution.


Case 3:
Lower Triangular Matrix :
   1.000    0.000    0.000 
   2.000    1.000    0.000 
   1.000    0.000    1.000 

Upper Triangular Matrix :
   1.000    1.000    1.000 
   0.000    0.000    0.000 
   0.000    0.000    2.000 

The system has infinite solutions.



```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/A.%20Solution%20of%20Linear%20Equations/LU%20Decomposition)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [LU Decomposition - Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
- [LU Decomposition Method - GeeksforGeeks](https://www.geeksforgeeks.org/l-u-decomposition-system-linear-equations/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="iterative-methods"></a>
### 4. Iterative Methods: Jacobi and Gauss-Seidel methods

<a id="jacobi"></a>
#### (i) Jacobi Iterative Method

**Theory: Simultaneous Displacement**
The Jacobi method is the simplest iterative technique. It works by isolating the variable $x_i$ in the $i$-th equation. The unique characteristic of Jacobi is that it uses values from the **previous** iteration to calculate **all** new values. No new information is used until the next full cycle.

Because each update is independent of the others within the same step, this method is famously easy to parallelize on modern multi-core processors.

**Algorithm**
1.  **Arrangement:** Rewrite the system so that $x_1$ is on the left of equation 1, $x_2$ on the left of equation 2, etc.
    - *Note: The system must be diagonally dominant for guaranteed convergence.*
2.  **Guess:** Start with an initial guess $x^{(0)}$ (often all zeros).
3.  **Iterate:** For each variable $x_i$, compute the new value using only the old values from step $k$.
4.  **Stop:** Repeat until the difference between the new and old values is less than your allowed error tolerance ($\epsilon$).

**Pseudocode**
```text
Input: Matrix A, Vector b, tolerance e, max_iterations N
Initialize: x_old = [0, 0, ... 0]
            x_new = [0, 0, ... 0]

For k from 1 to N:
    For i from 1 to rows(A):
        sum = 0
        For j from 1 to columns(A):
            if i != j:
                sum = sum + A[i][j] * x_old[j]
        
        x_new[i] = (b[i] - sum) / A[i][i]
    
    If distance(x_new, x_old) < e:
        Print "Converged"
        Break
    
    x_old = x_new // Update for next cycle
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream infile("inputjacobi.txt");
    ofstream outfile("outputjacobi.txt");

    if (!infile || !outfile)
        return 1;
    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);
        vector<double> x(n, 0.0);
        vector<double> new_x(n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                infile >> a[i][j];
            }
            infile >> b[i];
        }

        double tol;
        int max_iter;
        infile >> tol >> max_iter;

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Iteration\t";
        for (int i = 0; i < n; i++)
            outfile << "x" << i + 1 << "\t\t";
        outfile << "Error" << endl;
        outfile << endl;

        int iter = 0;
        double error = 0.0;

        do
        {
            error = 0.0;
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += a[i][j] * x[j];
                    }
                }
                new_x[i] = (b[i] - sum) / a[i][i];
                double diff = abs(new_x[i] - x[i]);
                if (diff > error)
                {
                    error = diff;
                }
            }

            x = new_x;
            iter++;

            outfile << iter << "\t\t";
            for (int i = 0; i < n; i++)
                outfile << fixed << setprecision(3) << x[i] << "\t";
            outfile << scientific << setprecision(3) << error << endl;

        } while (error > tol && iter < max_iter);

        outfile << endl;
        if (iter < max_iter)
        {
            outfile << "Converged in " << iter << " iterations." << endl;
        }
        else
        {
            outfile << "Maximum iterations reached." << endl;
        }

        outfile << "Solution:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << "x" << i + 1 << " = " << fixed << setprecision(3) << x[i] << endl;
        }
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
3
3 -0.1 -0.2 7.85
0.1 7 -0.3 -19.3
0.3 -0.2 10 71.4
0.00001
100
3
4 1 1 2
1 5 2 -6
1 2 4 -4
0.0001
100

```

#### Output
```
Case 1:
Iteration	x1		x2		x3		Error

1		2.617	-2.757	7.140	7.140e+000
2		3.001	-2.489	7.006	3.841e-001
3		3.001	-2.500	7.000	1.121e-002
4		3.000	-2.500	7.000	7.839e-004
5		3.000	-2.500	7.000	2.385e-005
6		3.000	-2.500	7.000	1.266e-006

Converged in 6 iterations.
Solution:
x1 = 3.000
x2 = -2.500
x3 = 7.000

Case 2:
Iteration	x1		x2		x3		Error

1		0.500	-1.200	-1.000	1.200e+000
2		1.050	-0.900	-0.525	5.500e-001
3		0.856	-1.200	-0.812	3.000e-001
4		1.003	-1.046	-0.614	1.984e-001
5		0.915	-1.155	-0.728	1.136e-001
6		0.971	-1.092	-0.651	7.639e-002
7		0.936	-1.134	-0.697	4.542e-002
8		0.958	-1.108	-0.667	2.955e-002
9		0.944	-1.125	-0.685	1.801e-002
10		0.952	-1.115	-0.674	1.151e-002
11		0.947	-1.121	-0.681	7.107e-003
12		0.950	-1.117	-0.676	4.496e-003
13		0.948	-1.120	-0.679	2.796e-003
14		0.950	-1.118	-0.677	1.760e-003
15		0.949	-1.119	-0.678	1.099e-003
16		0.949	-1.118	-0.678	6.899e-004
17		0.949	-1.119	-0.678	4.313e-004
18		0.949	-1.119	-0.678	2.705e-004
19		0.949	-1.119	-0.678	1.693e-004
20		0.949	-1.119	-0.678	1.061e-004
21		0.949	-1.119	-0.678	6.643e-005

Converged in 21 iterations.
Solution:
x1 = 0.949
x2 = -1.119
x3 = -0.678


```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/A.%20Solution%20of%20Linear%20Equations/Iterative%20Methods/Jacobi)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Jacobi Method - Wikipedia (Comprehensive Theory)](https://en.wikipedia.org/wiki/Jacobi_method)
* [Jacobi Method Explained - GeeksforGeeks (Examples & Code)](https://www.geeksforgeeks.org/engineering-mathematics/jacobian-method/)
* Numerical Methods for Engineers - Chapra & Canale

<a id="gaussseidel"></a>
#### (ii) Gauss-Seidel Iterative Method

**Theory: Successive Displacement**
The Gauss-Seidel method is an optimization of the Jacobi technique. In the Jacobi method, we hold all updates in a buffer until the end of the iteration. In Gauss-Seidel, we update the variables **immediately**.

As soon as a new value for $x_1$ is calculated, it replaces the old $x_1$ in memory and is used immediately to calculate $x_2$. Because we are using the most recent information available at every step, Gauss-Seidel typically converges much faster than Jacobi (often requiring half as many iterations). However, because $x_2$ depends on the *new* $x_1$, the steps must be done in order, making parallelization difficult.

**Algorithm**
1.  **Rearrange:** Same setup as Jacobi (requires diagonal dominance).
2.  **Initialize:** Start with an initial guess vector $x$.
3.  **Iterate:** Calculate $x_i$ and immediately overwrite the value in the solution vector. Use this new value for calculations of $x_{i+1}$, $x_{i+2}$, etc.
4.  **Check Convergence:** Stop when the error falls below the tolerance level.

**Pseudocode**
```text
Input: Matrix A, Vector b, tolerance e, max_iterations N
Output: Solution vector x

Initialize: 
    x = [0, 0, ... 0] // Only one array is needed

For k from 1 to N:
    converged = true
    
    For i from 1 to rows(A):
        old_val = x[i]
        sum = 0
        For j from 1 to columns(A):
            if i != j:
                // Uses the most recent x[j] available in memory
                sum = sum + A[i][j] * x[j] 
        
        x[i] = (b[i] - sum) / A[i][i]
        
        // Check error for this specific variable
        If abs(x[i] - old_val) > e:
            converged = false
            
    If converged is true:
        Print "Converged"
        Return x
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream infile("inputseidel.txt");
    ofstream outfile("outputseidel.txt");

    if (!infile)
        return 1;
    if (!outfile)
        return 1;

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);
        vector<double> x(n, 0.0);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                infile >> a[i][j];
            }
            infile >> b[i];
        }

        double tol;
        int max_iter;
        infile >> tol >> max_iter;

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Iteration\t";
        for (int i = 0; i < n; i++)
            outfile << "x" << i + 1 << "\t\t";
        outfile << "Error" << endl;
        outfile << endl;

        int iter = 0;
        double error = 0.0;

        do
        {
            error = 0.0;
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += a[i][j] * x[j];
                    }
                }
                double new_val = (b[i] - sum) / a[i][i];
                double diff = abs(new_val - x[i]);
                if (diff > error)
                {
                    error = diff;
                }
                x[i] = new_val;
            }

            iter++;

            outfile << iter << "\t\t";
            for (int i = 0; i < n; i++)
                outfile << fixed << setprecision(6) << x[i] << "\t";
            outfile << scientific << setprecision(4) << error << endl;

        } while (error > tol && iter < max_iter);

        outfile << endl;
        if (iter < max_iter)
        {
            outfile << "Converged in " << iter << " iterations." << endl;
        }
        else
        {
            outfile << "Maximum iterations reached." << endl;
        }

        outfile << "Solution:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << "x" << i + 1 << " = " << fixed << setprecision(3) << x[i] << endl;
        }
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
3
3 -0.1 -0.2 7.85
0.1 7 -0.3 -19.3
0.3 -0.2 10 71.4
0.00001
100
3
4 1 1 2
1 5 2 -6
1 2 4 -4
0.0001
100

```

#### Output
```
Case 1:
Iteration	x1		x2		x3		Error

1		2.616667	-2.794524	7.005610	7.0056e+000
2		2.990557	-2.499625	7.000291	3.7389e-001
3		3.000032	-2.499988	6.999999	9.4754e-003
4		3.000000	-2.500000	7.000000	3.1545e-005
5		3.000000	-2.500000	7.000000	3.5441e-007

Converged in 5 iterations.
Solution:
x1 = 3.000
x2 = -2.500
x3 = 7.000

Case 2:
Iteration	x1		x2		x3		Error

1		0.500000	-1.300000	-0.475000	1.3000e+000
2		0.943750	-1.198750	-0.636563	4.4375e-001
3		0.958828	-1.137141	-0.671137	6.1609e-002
4		0.952069	-1.121959	-0.677038	1.5181e-002
5		0.949749	-1.119135	-0.677870	2.8244e-003
6		0.949251	-1.118702	-0.677962	4.9806e-004
7		0.949166	-1.118649	-0.677967	8.5190e-005

Converged in 7 iterations.
Solution:
x1 = 0.949
x2 = -1.119
x3 = -0.678


```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/A.%20Solution%20of%20Linear%20Equations/Iterative%20Methods/GaussSeidel)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Gauss-Seidel Method - GeeksforGeeks](https://www.geeksforgeeks.org/gauss-seidel-method/)
* [Iterative Methods: Jacobi vs Gauss-Seidel - Math LibreTexts](https://math.libretexts.org/Bookshelves/Linear_Algebra/Introduction_to_Matrix_Algebra_(Kaw)/01%3A_Chapters/1.08%3A_Gauss-Seidel_Method)
* Numerical Methods for Engineers - Chapra & Canale



---

## B. Solution of Non-linear Equations

<a id="bisection"></a>
### 1. Bisection Method

**Theory**

The Bisection Method is a simple and robust numerical technique used to find roots of continuous functions. It is also known as the Binary Chopping Method or Half-Interval Method.

Given two points $x_1$ and $x_2$ where $f(x_1) \cdot f(x_2) < 0$, the midpoint approximation is:

$$x_0 = \frac{x_1 + x_2}{2}$$

The method repeatedly bisects the interval and selects the subinterval where the root lies.

**Algorithm**

1. **Choose interval:** Select $x_1$ and $x_2$ such that $f(x_1) \cdot f(x_2) < 0$
2. **Compute midpoint:** Calculate $x_0 = \frac{x_1 + x_2}{2}$
3. **Evaluate:** Compute $f(x_0)$
4. **Update interval:**
   - If $f(x_0) = 0$, then $x_0$ is the exact root
   - If $f(x_0) \cdot f(x_1) < 0$, set $x_2 = x_0$
   - If $f(x_0) \cdot f(x_2) < 0$, set $x_1 = x_0$
5. **Check convergence:** Repeat until $|x_2 - x_1| < \epsilon$
6. **Output:** The approximate root is $x_0$

**Pseudocode:**
```text
BisectionMethod(f, a, b, tolerance, max_iter):

  Step 1: Check if: f(a).f(b) < 0
          If not, print "Invalid interval" and stop.
  Step 2: For i = 1 to max_iter:
            x0 = (a + b) / 2
            fx0 = f(x0)
            If |fx0| < tolerance:
                Return x0 as the root
            Else If f(a) . fx0 < 0:
                b = x0
            Else:
                a = x0
            End For

  Step 3: Return (a + b) / 2 as the approximate root
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(const vector<double> &coeffs, double x)
{
    double result = 0.0;
    int n = coeffs.size();
    for (int i = 0; i < n; i++)
    {
        result += coeffs[i] * pow(x, n - i - 1);
    }
    return result;
}

bool stop(double prev, double curr, double tol)
{
    return fabs(prev - curr) < tol;
}

vector<double> Bisection(const vector<double> &coeffs, double l, double u, double tol, double step)
{
    vector<pair<double, double>> intervals;
    vector<double> roots;

    for (double x = l; x <= u - step; x += step)
    {
        double fx = f(coeffs, x);
        double fx_next = f(coeffs, x + step);

        if (fabs(fx) < 1e-12)
        {
            roots.push_back(x);
            continue;
        }

        if (fx * fx_next < 0)
        {
            intervals.push_back({x, x + step});
        }
    }

    for (auto interval : intervals)
    {
        double a = interval.first;
        double b = interval.second;
        double fa = f(coeffs, a);
        double fb = f(coeffs, b);

        if (fa * fb > 0)
            continue;

        double prev = a;
        double mid = (a + b) / 2.0;
        int iter = 0;

        while (!stop(prev, mid, tol))
        {
            double fmid = f(coeffs, mid);
            if (fabs(fmid) < 1e-12)
                break;

            if (fa * fmid < 0)
            {
                b = mid;
                fb = fmid;
            }
            else
            {
                a = mid;
                fa = fmid;
            }

            prev = mid;
            mid = (a + b) / 2.0;
            iter++;
        }

        roots.push_back(mid);
    }

    sort(roots.begin(), roots.end());
    roots.erase(unique(roots.begin(), roots.end(), [tol](double a, double b)
                       { return fabs(a - b) < tol; }),
                roots.end());

    return roots;
}

string print_function(const vector<double> &coeffs)
{
    stringstream out;
    int n = coeffs.size();
    out << "f(x)=";
    for (int i = 0; i < n; i++)
    {
        if (i != 0 && coeffs[i] >= 0)
            out << "+";
        out << coeffs[i] << "x^" << n - i - 1 << " ";
    }
    return out.str();
}

int main()
{
    ifstream infile("inputbisection.txt");
    ofstream outfile("outputbisection.txt");

    if (!infile)
    {
        return 1;
    }

    int numProblems;
    infile >> numProblems;

    for (int p = 0; p < numProblems; p++)
    {
        int degree;
        infile >> degree;

        vector<double> coeffs(degree + 1);
        for (int i = 0; i <= degree; i++)
        {
            infile >> coeffs[i];
        }

        double l, u, tol, step;
        infile >> l >> u >> tol >> step;

        outfile << "Problem " << p + 1 << endl;
        outfile << "Function: " << print_function(coeffs) << endl;
        outfile << "Interval: [" << l << ", " << u << "], Tolerance: " << tol << ", Step: " << step << endl;

        vector<double> roots = Bisection(coeffs, l, u, tol, step);

        outfile << "Roots found : " << endl;
        if (!roots.empty())
        {
            for (double root : roots)
            {
                outfile << fixed << setprecision(6) << root << endl;
            }
        }
        outfile << "\n\n"
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
3
2
1 -4 -10
5 6 0.0001 0.01
3
1 0 -1 -2
-2 3 0.0001 0.01
5
1 -15 85 -225 274 -120
0 5 0.0001 0.01
```

#### Output
```
Problem 1
Function: f(x)=1x^2 -4x^1 -10x^0 
Interval: [5, 6], Tolerance: 0.0001, Step: 0.01
Roots found : 
5.741641



Problem 2
Function: f(x)=1x^3 +0x^2 -1x^1 -2x^0 
Interval: [-2.000000, 3.000000], Tolerance: 0.000100, Step: 0.010000
Roots found : 
1.521328



Problem 3
Function: f(x)=1x^5 -15x^4 +85x^3 -225x^2 +274x^1 -120x^0 
Interval: [0.000000, 5.000000], Tolerance: 0.000100, Step: 0.010000
Roots found : 
0.999922
1.999922
3.000000
4.000000




```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/B.%20Solution%20of%20Non-Linear%20Equations/Bisection)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Bisection Method - Wikipedia](https://en.wikipedia.org/wiki/Bisection_method)
- Numerical Methods for Engineers - Chapra & Canale

<a id="false_position"></a>
### 2. False Position Method

**Theory**

The False Position Method (also known as Regula Falsi) is a bracketing method similar to Bisection, but instead of using the midpoint, it uses a linear interpolation to find a better approximation of the root.

Given two points $x_1$ and $x_2$ where $f(x_1) \cdot f(x_2) < 0$, the next approximation is:

$$x_0 = x_1 - \frac{f(x_1) \cdot (x_2 - x_1)}{f(x_2) - f(x_1)}$$

This method converges faster than Bisection because it uses the function values to estimate where the root lies, rather than simply bisecting the interval.

**Algorithm of False Position Method:**

Step 1: Choosing two numbers x1 and x2 such that:
    f(x1).f(x2)<0
Step 2: Compute the root approximation:
    x0 = x1 - (f(x1).(x2-x1))/(f(x2)-f(x1))
Step 3: Evaluate f(x0).
Step 4:
If f(x0)=0, then x0 is the exact root.
If f(x0).f(x1)<0  then, x2=x0
If f(x0).f(x2)<0  then, x1=x0
Step 5: Repeat until the stopping criterion is met:
|x2-x1|<E
Step 6: Stop. The approximate root is x0

**Pseudocode:**
```text
FalsePositionMethod(f, a, b, tolerance, max_iter):

  Step 1: Check if f(a) . f(b) < 0
          If not, print "Invalid interval" and stop.
  Step 2: For i = 1 to max_iter:
            x0 = a - f(a) * (b - a) / (f(b) - f(a))
            fx0 = f(x0)
            If |fx0| < tolerance:
                Return x0 as the root
            Else If f(a) . fx0 < 0:
                b = x0
            Else:
                a = x0
  Step 3: Return x0 as the approximate root
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(const vector<double> &coeffs, double x)
{
    double result = 0.0;
    int n = coeffs.size();
    for (int i = 0; i < n; i++)
    {
        result += coeffs[i] * pow(x, n - i - 1);
    }
    return result;
}

bool stop(double prev, double curr, double tol)
{
    return fabs(prev - curr) < tol;
}

vector<double> FalsePosition(const vector<double> &coeffs, double l, double u, double tol, double step)
{
    vector<pair<double, double>> intervals;
    vector<double> roots;

    for (double x = l; x <= u - step; x += step)
    {
        double fx = f(coeffs, x);
        double fx_next = f(coeffs, x + step);

        if (fabs(fx) < 1e-12)
        {
            roots.push_back(x);
            continue;
        }

        if (fx * fx_next < 0)
        {
            intervals.push_back({x, x + step});
        }
    }

    for (auto interval : intervals)
    {
        double a = interval.first;
        double b = interval.second;
        double fa = f(coeffs, a);
        double fb = f(coeffs, b);

        if (fa * fb > 0)
            continue;

        double prev = a;
        double c = a - (fa * (b - a)) / (fb - fa);
        int iter = 0;

        while (!stop(prev, c, tol))
        {
            double fc = f(coeffs, c);
            if (fabs(fc) < 1e-12)
                break;

            if (fa * fc < 0)
            {
                b = c;
                fb = fc;
            }
            else
            {
                a = c;
                fa = fc;
            }

            prev = c;
            c = a - (fa * (b - a)) / (fb - fa);
            iter++;
        }

        roots.push_back(c);
    }

    sort(roots.begin(), roots.end());
    roots.erase(unique(roots.begin(), roots.end(), [tol](double a, double b)
                       { return fabs(a - b) < tol; }),
                roots.end());

    return roots;
}

string print_function(const vector<double> &coeffs)
{
    stringstream out;
    int n = coeffs.size();
    out << "f(x)=";
    for (int i = 0; i < n; i++)
    {
        if (i != 0 && coeffs[i] >= 0)
            out << "+";
        out << coeffs[i] << "x^" << n - i - 1 << " ";
    }
    return out.str();
}

int main()
{
    ifstream infile("inputfalse.txt");
    ofstream outfile("outputfalse.txt");

    if (!infile)
    {
        return 1;
    }

    int numProblems;
    infile >> numProblems;

    for (int p = 0; p < numProblems; p++)
    {
        int degree;
        infile >> degree;

        vector<double> coeffs(degree + 1);
        for (int i = 0; i <= degree; i++)
        {
            infile >> coeffs[i];
        }

        double l, u, tol, step;
        infile >> l >> u >> tol >> step;

        outfile << "Problem " << p + 1 << endl;
        outfile << "Function: " << print_function(coeffs) << endl;
        outfile << "Interval: [" << l << ", " << u << "], Tolerance: " << tol << ", Step: " << step << endl;

        vector<double> roots = FalsePosition(coeffs, l, u, tol, step);

        outfile << "Roots found : " << endl;
        if (!roots.empty())
        {
            for (double root : roots)
            {
                outfile << fixed << setprecision(6) << root << endl;
            }
        }
        outfile << "\n\n"
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
1
5 1 -15 85 -225 274 -120 0 6 0.0001 0.01
```

#### Output
```
Problem 1
Function: f(x)=1x^5 -15x^4 +85x^3 -225x^2 +274x^1 -120x^0 
Interval: [0, 6], Tolerance: 0.0001, Step: 0.01
Roots found : 
1.000000
2.000000
3.000000
4.000000
5.000000

```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/B.%20Solution%20of%20Non-Linear%20Equations/False_Position)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [False Position Method - Wikipedia](https://en.wikipedia.org/wiki/Regula_falsi)
- [False Position Method - GeeksforGeeks](https://www.geeksforgeeks.org/program-for-method-of-false-position/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="secant"></a>
### 3. Secant Method

**Theory**

The Secant Method approximates the derivative numerically using two points instead of requiring an analytical derivative. Given two initial points $x_{n-1}$ and $x_n$:

$$f'(x_n) \approx \frac{f(x_n) - f(x_{n-1})}{x_n - x_{n-1}}$$

Substituting this into Newton's formula gives the Secant iteration:

$$x_{n+1} = x_n - \frac{f(x_n) \cdot (x_n - x_{n-1})}{f(x_n) - f(x_{n-1})}$$

**Algorithm**

1. **Initialize:** Choose two initial approximations $x_0$ and $x_1$
2. **Iterate:** Apply the Secant formula to compute $x_{n+1}$
3. **Check convergence:** If $|x_{n+1} - x_n| < \epsilon$, stop
4. **Update:** Set $x_{n-1} = x_n$ and $x_n = x_{n+1}$, then repeat

**Pseudocode**
```text
SecantMethod(f, x0, x1, tolerance, max_iter):
for i = 1 to max_iter:
    f0 = f(x0)
    f1 = f(x1)
    x2 = x1 - f1*(x1 - x0)/(f1 - f0)
    if abs(x2 - x1) < tolerance:
        return x2
    x0 = x1
    x1 = x2
return x2
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> coef;
int degree;

string polynomial_string(const vector<double> &c, int deg)
{
    stringstream ss;
    bool first = true;
    for (int i = 0; i <= deg; i++)
    {
        if (c[i] == 0)
            continue;
        if (!first && c[i] > 0)
            ss << "+";
        ss << c[i];
        if (deg - i > 0)
            ss << "x";
        if (deg - i > 1)
            ss << "^" << deg - i;
        first = false;
    }
    return ss.str();
}

double f(double x)
{
    double value = 0.0;
    for (int i = 0; i <= degree; i++)
        value += coef[i] * pow(x, degree - i);
    return value;
}

double df(double x)
{
    double value = 0.0;
    for (int i = 0; i < degree; i++)
        value += coef[i] * (degree - i) * pow(x, degree - i - 1);
    return value;
}

vector<pair<double, double>> initial_guess_pairs(double step)
{
    double dmax = 0.0;
    for (int j = 0; j <= degree; j++)
        dmax = max(dmax, fabs(coef[j] / coef[0]));

    double low = -(1 + dmax);
    double high = (1 + dmax);
    vector<pair<double, double>> res;

    while (low < high)
    {
        double x1 = low, x2 = low + step;
        if (f(x1) * f(x2) < 0)
            res.push_back({x1, x2});
        low += step;
    }
    return res;
}

double secant_method(double x0, double x1, double tol)
{
    double f0 = f(x0);
    double f1 = f(x1);

    int iter = 0;
    while (fabs(x1 - x0) >= tol)
    {
        if (fabs(f1 - f0) < 1e-12)
            return NAN;
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);

        iter++;
        if (iter > 1000)
            break;
    }
    return x1;
}

int main()
{
    ifstream fin("inputsecant.txt");
    ofstream fout("outputsecant.txt");

    int T;
    fin >> T; 

    for (int t = 1; t <= T; t++)
    {
        fin >> degree;
        coef.resize(degree + 1);
        for (int i = 0; i <= degree; i++)
            fin >> coef[i];

        fout << "Polynomial " << t << ":\n";
        fout << "f(x) = " << polynomial_string(coef, degree) << "\n";

 
        vector<double> dcoef(degree);
        for (int i = 0; i < degree; i++)
            dcoef[i] = coef[i] * (degree - i);
        fout << "f'(x) = " << polynomial_string(dcoef, degree - 1) << "\n";

        vector<pair<double, double>> guesses = initial_guess_pairs(0.45);
        double tol = 1e-6;
        vector<double> roots;

        for (auto &p : guesses)
        {
            double root = secant_method(p.first, p.second, tol);
            if (isnan(root))
                continue;

            bool unique = true;
            for (double r : roots)
                if (fabs(root - r) < tol)
                    unique = false;

            if (unique)
                roots.push_back(root);
        }

        fout << "Roots:\n";
        for (double r : roots)
            fout << r << "\n";
        fout << "\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Input
```
2
3
1 -6 11 -6
2
1 -3 2
```

#### Output
```
Polynomial 1:
f(x) = 1x^3-6x^2+11x-6
f'(x) = 3x^2-12x+11
Roots:
1
2
3


Polynomial 2:
f(x) = 1x^2-3x+2
f'(x) = 2x-3
Roots:
1
2



```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/B.%20Solution%20of%20Non-Linear%20Equations/Secant)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Secant Method - Wikipedia](https://en.wikipedia.org/wiki/Secant_method)
- [Secant Method - GeeksforGeeks](https://www.geeksforgeeks.org/program-to-find-root-of-an-equations-using-secant-method/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="newton_raphson"></a>
### 4. Newton Raphson Method

**Theory**

The Newton-Raphson method uses the tangent line at the current approximation to estimate a better root approximation. Given a guess $x_n$, the tangent line at that point is:

$$y = f(x_n) + f'(x_n)(x - x_n)$$

Setting $y = 0$ and solving for $x$:

$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$

This is the Newton-Raphson formula.

**Algorithm**

1. **Initialize:** Choose an initial guess $x_0$
2. **Compute derivative:** Calculate $f'(x_n)$
3. **Compute next approximation:** Apply $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
4. **Check convergence:** If $|x_{n+1} - x_n| < \epsilon$, stop
5. **Update:** Set $x_n = x_{n+1}$ and repeat

**Pseudocode**
```text
NewtonRaphson(f, df, x0, tolerance, max_iter):
for i = 1 to max_iter:
    fx = f(x0)
    dfx = df(x0)
    if dfx == 0:
        break
    x1 = x0 - fx/dfx
    if abs(x1 - x0) < tolerance:
        return x1
    x0 = x1
return x1
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> coef;
int degree;

string polynomial_string(const vector<double> &c, int deg)
{
    stringstream ss;
    bool first = true;
    for (int i = 0; i <= deg; i++)
    {
        if (c[i] == 0)
            continue;
        if (!first && c[i] > 0)
            ss << "+";
        ss << c[i];
        if (deg - i > 0)
            ss << "x";
        if (deg - i > 1)
            ss << "^" << deg - i;
        first = false;
    }
    return ss.str();
}

double f(double x)
{
    double value = 0.0;
    for (int i = 0; i <= degree; i++)
        value += coef[i] * pow(x, degree - i);
    return value;
}

double df(double x)
{
    double value = 0.0;
    for (int i = 0; i < degree; i++)
        value += coef[i] * (degree - i) * pow(x, degree - i - 1);
    return value;
}

vector<double> rangess(double step)
{
    double dmax = 0.0;
    for (int j = 0; j <= degree; j++)
        dmax = max(dmax, fabs(coef[j] / coef[0]));

    double low = -(1 + dmax);
    double high = (1 + dmax);
    vector<double> res;

    while (low < high)
    {
        if (f(low) * f(low + step) < 0)
        {
            if (fabs(f(low)) < fabs(f(low + step)))
                res.push_back(low);
            else
                res.push_back(low + step);
        }
        low += step;
    }
    return res;
}

double newton_raphson(double initial_guess, double tol)
{
    double x0 = initial_guess;
    double fx = f(x0);
    if (fabs(fx) < 1e-12)
        return x0;

    double dfx = df(x0);
    if (fabs(dfx) <= 1e-12)
        return NAN;

    double x1 = x0 - fx / dfx;

    while (fabs(x1 - x0) >= tol)
    {
        x0 = x1;
        fx = f(x0);
        if (fabs(fx) < 1e-12)
            return x0;
        dfx = df(x0);
        if (fabs(dfx) <= 1e-12)
            return NAN;
        x1 = x0 - fx / dfx;
    }
    return x1;
}

int main()
{
    ifstream fin("inputnewtraph.txt");
    ofstream fout("outputnewtraph.txt");

    int T;
    fin >> T;

    for (int t = 1; t <= T; t++)
    {
        fin >> degree;
        coef.resize(degree + 1);
        for (int i = 0; i <= degree; i++)
            fin >> coef[i];

        fout << "Polynomial " << t << ":\n";
        fout << "f(x) = " << polynomial_string(coef, degree) << "\n";

        vector<double> dcoef(degree);
        for (int i = 0; i < degree; i++)
            dcoef[i] = coef[i] * (degree - i);
        fout << "f'(x) = " << polynomial_string(dcoef, degree - 1) << "\n";

        vector<double> initial_guesses = rangess(0.45);
        double tol = 1e-6;
        vector<double> roots;

        for (double guess : initial_guesses)
        {
            double root = newton_raphson(guess, tol);
            if (isnan(root))
                continue;

            bool unique = true;
            for (double r : roots)
                if (fabs(root - r) < tol)
                    unique = false;

            if (unique)
                roots.push_back(root);
        }

        fout << "Roots:\n";
        for (double r : roots)
            fout << r << "\n";

        fout << "\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Input
```
2
3
1 -6 11 -6
2
1 -3 2
```

#### Output
```
Polynomial 1:
f(x) = 1x^3-6x^2+11x-6
f'(x) = 3x^2-12x+11
Roots:
1
2
3


Polynomial 2:
f(x) = 1x^2-3x+2
f'(x) = 2x-3
Roots:
1
2


```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/B.%20Solution%20of%20Non-Linear%20Equations/Newton_Raphson)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Newton-Raphson Method - Wikipedia](https://en.wikipedia.org/wiki/Newton%27s_method)
- [Newton-Raphson Method - GeeksforGeeks](https://www.geeksforgeeks.org/newton-raphson-method/)
- Numerical Methods for Engineers - Chapra & Canale

---

## C. Interpolation and Approximation

<a id="newton-forward"></a>
### 1. Newton Forward Interpolation Method

**Theory**
Newton's Forward Interpolation is used to approximate the value of a function $f(x)$ at valid points, based on a set of known data points that are **equally spaced**.

The interpolation formula is given by:
$$y(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \frac{u(u-1)(u-2)}{3!} \Delta^3 y_0 + \dots$$
Where $u = \frac{x - x_0}{h}$.

**Algorithm**
1.  **Input**: Read ($n$) data points ($x, y$) and the value to interpolate ($value$).
2.  **Difference Table**: Construct the forward difference table.
3.  **Calculate u**: Compute $u = (value - x[0]) / (x[1] - x[0])$.
4.  **Compute Sum**: Apply formula.
5.  **Output**: Display the interpolated value.

**Pseudocode**
```cpp
function newton_forward(x[], y[][], n, value):
    for i from 1 to n-1:
        for j from 0 to n-i-1:
            y[j][i] = y[j+1][i-1] - y[j][i-1]

    sum = y[0][0]
    u = (value - x[0]) / (x[1] - x[0])
    
    for i from 1 to n-1:
        u_term = 1
        for k from 0 to i-1:
            u_term = u_term * (u - k) 
        
        factorial = fact(i)
        sum = sum + (u_term * y[0][i]) / factorial

    return sum
```
#### Code
```cpp
#include <bits/stdc++.h>

using namespace std;

double calculateU(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u + i);
    return temp;
}

double fact(int n)
{
    double f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream infile("inputbackward.txt");
    ofstream outfile("outputbackward.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n);
        vector<vector<double>> y(n, vector<double>(n));

        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i][0];
        }

        for (int i = 1; i < n; i++)
        {
            for (int j = n - 1; j >= i; j--)
                y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Backward Difference Table:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << setw(10) << x[i];
            for (int j = 0; j <= i; j++)
                outfile << setw(10) << fixed << setprecision(4) << y[i][j];
            outfile << endl;
        }

        double value;
        infile >> value;

        double sum = y[n - 1][0];
        double u = (value - x[n - 1]) / (x[1] - x[0]);
        double last_term = 0.0;

        for (int i = 1; i < n; i++)
        {
            last_term = (calculateU(u, i) * y[n - 1][i]) / fact(i);
            sum = sum + last_term;
        }

        outfile << endl
                << "Value at " << value << " is " << fixed << setprecision(6) << sum << endl;

        if (sum != 0)
        {
            double relative_error = abs(last_term / sum) * 100.0;
            outfile << "Approximate Relative Error: " << relative_error << "%" << endl;
        }
        else
        {
            outfile << "Approximate Error: " << abs(last_term) << " (Cannot calculate percentage relative error as sum is 0)" << endl;
        }

        outfile << endl
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
4
45 0.7071
50 0.7660
55 0.8192
60 0.8660
46
5
1891 46
1901 66
1911 81
1921 93
1931 101
1895


```

#### Output
```
Case 1:
Difference Table:
        45    0.7071    0.0589   -0.0057   -0.0007
   50.0000    0.7660    0.0532   -0.0064
   55.0000    0.8192    0.0468
   60.0000    0.8660

Value at 46.0000 is 0.719

Case 2:
Difference Table:
  1891.000   46.0000   20.0000   -5.0000    2.0000   -3.0000
 1901.0000   66.0000   15.0000   -3.0000   -1.0000
 1911.0000   81.0000   12.0000   -4.0000
 1921.0000   93.0000    8.0000
 1931.0000  101.0000

Value at 1895.0000 is 54.853





```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/C.%20Interpolation%20and%20Approximation/Newton%20Forward)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Newton Forward Interpolation - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- Numerical Methods for Engineers - Chapra & Canale

<a id="newton-backward"></a>
### 2. Newton Backward Interpolation Method

**Theory**
Newton's Backward Interpolation is similar to the forward method but is more accurate for interpolating values near the **end** of the dataset.

The formula is given by:
$$y(x) = y_n + u \nabla y_n + \frac{u(u+1)}{2!} \nabla^2 y_n + \frac{u(u+1)(u+2)}{3!} \nabla^3 y_n + \dots$$
Where $u = \frac{x - x_n}{h}$.

**Algorithm**
1.  **Input**: Read ($n$) data points ($x, y$) and the values.
2.  **Difference Table**: Construct the backward difference table.
3.  **Calculate u**: Compute $u = (value - x[n-1]) / (x[1] - x[0])$.
4.  **Compute Sum**: Apply formula.
5.  **Output**: Display result.

**Pseudocode**
```cpp
function newton_backward(x[], y[][], n, value):
    for i from 1 to n-1:
        for j from n-1 down to i:
             y[j][i] = y[j][i-1] - y[j-1][i-1]

    sum = y[n-1][0]
    u = (value - x[n-1]) / (x[1] - x[0])

    for i from 1 to n-1:
        u_term = 1
        for k from 0 to i-1:
             u_term = u_term * (u + k) 
             
        factorial = fact(i)
        sum = sum + (u_term * y[n-1][i]) / factorial

    return sum
```
#### Code
```cpp
#include <bits/stdc++.h>

using namespace std;

double calculateU(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u + i);
    return temp;
}

double fact(int n)
{
    double f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream infile("inputbackward.txt");
    ofstream outfile("outputbackward.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n);
        vector<vector<double>> y(n, vector<double>(n));

        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i][0];
        }

        for (int i = 1; i < n; i++)
        {
            for (int j = n - 1; j >= i; j--)
                y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Backward Difference Table:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << setw(10) << x[i];
            for (int j = 0; j <= i; j++)
                outfile << setw(10) << fixed << setprecision(4) << y[i][j];
            outfile << endl;
        }

        double value;
        infile >> value;

        double sum = y[n - 1][0];
        double u = (value - x[n - 1]) / (x[1] - x[0]);
        double last_term = 0.0;

        for (int i = 1; i < n; i++)
        {
            last_term = (calculateU(u, i) * y[n - 1][i]) / fact(i);
            sum = sum + last_term;
        }

        outfile << endl
                << "Value at " << value << " is " << fixed << setprecision(6) << sum << endl;

        if (sum != 0)
        {
            double relative_error = abs(last_term / sum) * 100.0;
            outfile << "Approximate Relative Error: " << relative_error << "%" << endl;
        }
        else
        {
            outfile << "Approximate Error: " << abs(last_term) << " (Cannot calculate percentage relative error as sum is 0)" << endl;
        }

        outfile << endl
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
5
35 31
45 73
55 124
65 159
75 190
70
4
0 1
1 2
2 1
3 10
2.5

```

#### Output
```
Case 1:
Backward Difference Table:
        35   31.0000
   45.0000   73.0000   42.0000
   55.0000  124.0000   51.0000    9.0000
   65.0000  159.0000   35.0000  -16.0000  -25.0000
   75.0000  190.0000   31.0000   -4.0000   12.0000   37.0000

Value at 70.0000 is 172.804688
Approximate Relative Error: 0.836385%


Case 2:
Backward Difference Table:
  0.000000    1.0000
    1.0000    2.0000    1.0000
    2.0000    1.0000   -1.0000   -2.0000
    3.0000   10.0000    9.0000   10.0000   12.0000

Value at 2.5000 is 3.500000
Approximate Relative Error: 21.428571%



```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/C.%20Interpolation%20and%20Approximation/Newton%20Backward)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Newton Backward Interpolation - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- Numerical Methods for Engineers - Chapra & Canale

---


## D. Numerical Differentiation

<a id="newton-forward-differentiation"></a>
### 1. Newton Forward Differentiation

**Theory**

Newton's Forward Differentiation formula is used to compute the derivative of a function at a given point using a forward difference table. It is based on Newton's Forward Interpolation formula.

For equally spaced data points, the first derivative is approximated by:

$$f'(x) = \frac{1}{h} \left[ \Delta y_0 - \frac{2u-1}{2!} \Delta^2 y_0 + \frac{3u^2 - 6u + 2}{3!} \Delta^3 y_0 - \cdots \right]$$

Where:
- $h$ is the step size (interval between consecutive x values)
- $u = \frac{x - x_0}{h}$
- $\Delta y_0, \Delta^2 y_0, \ldots$ are forward differences

**Algorithm**

1. **Input:** Read $n$ data points $(x_i, y_i)$, the value at which derivative is needed, and the order of derivative.
2. **Difference Table:** Construct the forward difference table.
3. **Calculate $u$:** Compute $u = \frac{x - x_0}{h}$ where $h = x_1 - x_0$.
4. **Apply Formula:** Use the appropriate Newton's forward differentiation formula based on the order.
5. **Output:** Display the derivative value.

**Pseudocode**
```text
START

Read n
Read x[i], y[i][0]

FOR i = 1 to n-1
    FOR j = 0 to n-i-1
        y[j][i] = y[j+1][i-1] - y[j][i-1]
    END FOR
END FOR

Read value, order

h = x[1] - x[0]
u = (value - x[0]) / h

derivative = 0
P = 1

FOR k = 0 to n-1
    IF k > 0
        P = P × (u - (k-1))
    END IF

    IF k ≥ order
        Differentiate P (order times)
        derivative += (y[0][k] / k!) × P(u)
    END IF
END FOR

derivative = derivative / h^order

Print derivative

END
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;
typedef vector<double> Poly;

Poly multiply(const Poly &p, double constant)
{
    Poly res(p.size() + 1, 0.0);
    for (size_t i = 0; i < p.size(); i++)
    {
        res[i + 1] += p[i];
        res[i] += constant * p[i];
    }
    return res;
}

Poly differentiate(const Poly &p)
{
    if (p.size() <= 1)
        return {0.0};
    Poly res;
    for (size_t i = 1; i < p.size(); i++)
    {
        res.push_back(p[i] * i);
    }
    return res;
}

double evaluate(const Poly &p, double x)
{
    double res = 0;
    double x_pow = 1;
    for (double c : p)
    {
        res += c * x_pow;
        x_pow *= x;
    }
    return res;
}

double fact(int n)
{
    double f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream infile("inputdiff.txt");
    ofstream outfile("outputdiff.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n);
        vector<vector<double>> y(n, vector<double>(n));

        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i][0];
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = 0; j < n - i; j++)
                y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Forward Difference Table:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << setw(10) << x[i];
            for (int j = 0; j < n - i; j++)
                outfile << setw(10) << fixed << setprecision(4) << y[i][j];
            outfile << endl;
        }

        double value;
        int order;
        infile >> value >> order;

        double h = x[1] - x[0];
        double u = (value - x[0]) / h;

        double derivative_val = 0;
        Poly P = {1.0};

        for (int k = 0; k < n; k++)
        {
            if (k > 0)
            {
                P = multiply(P, -(double)(k - 1));
            }

            if (k >= order)
            {
                Poly P_deriv = P;
                for (int m = 0; m < order; m++)
                {
                    P_deriv = differentiate(P_deriv);
                }

                double term_val = evaluate(P_deriv, u);
                derivative_val += (y[0][k] / fact(k)) * term_val;
            }
        }

        derivative_val /= pow(h, order);

        outfile << endl
                << "For x = " << value << ":" << endl;
        outfile << "Derivative Order: " << order << endl;
        outfile << "Result: " << fixed << setprecision(3) << derivative_val << endl;
        outfile << endl
                << endl;
    }
    infile.close();
    outfile.close();

    return 0;
}

```

#### Input
```
7
1.0 2.7183
1.2 3.3201
1.4 4.0552
1.6 4.9530
1.8 6.0496
2.0 7.3891
2.2 9.0250
1.2 1
5
0 0
1 1
2 8
3 27
4 64
2 2

```

#### Output
```
Case 1:
Forward Difference Table:
         1    2.7183    0.6018    0.1333    0.0294    0.0067    0.0013    0.0001
    1.2000    3.3201    0.7351    0.1627    0.0361    0.0080    0.0014
    1.4000    4.0552    0.8978    0.1988    0.0441    0.0094
    1.6000    4.9530    1.0966    0.2429    0.0535
    1.8000    6.0496    1.3395    0.2964
    2.0000    7.3891    1.6359
    2.2000    9.0250

For x = 1.2000:
Derivative Order: 1
Result: 3.320


Case 2:
Forward Difference Table:
     0.000    0.0000    1.0000    6.0000    6.0000    0.0000
    1.0000    1.0000    7.0000   12.0000    6.0000
    2.0000    8.0000   19.0000   18.0000
    3.0000   27.0000   37.0000
    4.0000   64.0000

For x = 2.0000:
Derivative Order: 2
Result: 12.000





```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/E.%20Solution%20of%20Differential%20Equations/Newton%20Forward%20Differentiation)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Numerical Differentiation - Wikipedia](https://en.wikipedia.org/wiki/Numerical_differentiation)
- [Newton's Forward Difference Formula - MathWorld](https://mathworld.wolfram.com/ForwardDifference.html)
- Numerical Methods for Engineers - Chapra & Canale

  
## E. Solution of Differential Equations

<a id="runge-kutta"></a>
### 2. Runge Kutta Method

**Theory**
The Runge-Kutta method (specifically the fourth-order RK4) is a widely used technique for the approximate solution of ordinary differential equations (ODEs).
$$y_{n+1} = y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

**Algorithm**
1.  **Define function**: $f(x, y)$.
2.  **Input**: Initial $x_0, y_0$, target $x_n$, and step size $h$.
3.  **Iterate**: Until current $x$ reaches target $x_n$, calculate slopes $K$ and update $y$.
4.  **Output**: Final value of $y$.

**Pseudocode**
```cpp
function solve_rk4(x0, y0, xn, h):
    n_steps = (xn - x0) / h
    y = y0
    x = x0

    for i from 0 to n_steps:
        if x >= xn: break

        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y + k1/2)
        k3 = h * f(x + h/2, y + k2/2)
        k4 = h * f(x + h, y + k3)

        k_avg = (k1 + 2*k2 + 2*k3 + k4) / 6.0
        
        y = y + k_avg
        x = x + h
        
    return y
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x, double y)
{
    return x + y;
}

int main()
{
    ifstream infile("inputrk.txt");
    ofstream outfile("outputrk.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    double x0, y0, xn, h;
    int caseNum = 1;

    outfile << fixed << setprecision(3);

    while (infile >> x0 >> y0 >> xn >> h)
    {
        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Initial values: x0 = " << x0 << ", y0 = " << y0 << endl;
        outfile << "Target x: " << xn << endl;
        outfile << "Step size: " << h << endl;
        outfile << "List of x and y values:" << endl;
        outfile << "x\t\ty" << endl;

        double y = y0;
        double x = x0;
        outfile << x << "\t\t" << y << endl;
        int n = (int)round((xn - x0) / h);

        for (int i = 0; i < n; i++)
        {
            double k1 = h * f(x, y);
            double k2 = h * f(x + h / 2.0, y + k1 / 2.0);
            double k3 = h * f(x + h / 2.0, y + k2 / 2.0);
            double k4 = h * f(x + h, y + k3);

            double k = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

            y = y + k;
            x = x + h;

            outfile << x << "\t\t" << y << endl;
        }

        outfile << "Value of y at x = " << xn << " is " << y << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
0 1 0.2 0.1
0 1 0.5 0.1
1 2 1.5 0.1
```

#### Output
```
Case 1:
Initial values: x0 = 0.000, y0 = 1.000
Target x: 0.200
Step size: 0.100
List of x and y values:
x		y
0.000		1.000
0.100		1.110
0.200		1.243
Value of y at x = 0.200 is 1.243

Case 2:
Initial values: x0 = 0.000, y0 = 1.000
Target x: 0.500
Step size: 0.100
List of x and y values:
x		y
0.000		1.000
0.100		1.110
0.200		1.243
0.300		1.400
0.400		1.584
0.500		1.797
Value of y at x = 0.500 is 1.797

Case 3:
Initial values: x0 = 1.000, y0 = 2.000
Target x: 1.500
Step size: 0.100
List of x and y values:
x		y
1.000		2.000
1.100		2.321
1.200		2.686
1.300		3.099
1.400		3.567
1.500		4.095
Value of y at x = 1.500 is 4.095


```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/E.%20Solution%20of%20Differential%20Equations/Runge%20Kutta)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Runge-Kutta Methods - Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
- [Runge-Kutta Method - GeeksforGeeks](https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/)
- Numerical Methods for Engineers - Chapra & Canale

---

## F. Numerical Integration

<a id="simpson-1-by-3"></a>
### 1. Simpson’s 1/3 Rule

**Theory: Parabolic Approximation**
Simpson’s 1/3 Rule improves upon the Trapezoidal Rule by approximating the function $f(x)$ not as a straight line, but as a **second-order polynomial (parabola)** connecting every three points.



Because it fits parabolas, it requires an **even number of segments** (intervals) $n$, which means you need an odd number of data points. The formula weights the boundary points and internal points differently to achieve higher accuracy.

$$I \approx \frac{h}{3} \left[ (y_0 + y_n) + 4(y_1 + y_3 + \dots + y_{n-1}) + 2(y_2 + y_4 + \dots + y_{n-2}) \right]$$

**Algorithm**
1.  **Verify Intervals:** Ensure the number of intervals $n$ is even. If $n$ is odd, this method cannot be applied directly over the whole range.
2.  **Calculate Step Size:** $h = (b - a) / n$.
3.  **Sum Extremes:** Add the first ($y_0$) and last ($y_n$) values.
4.  **Sum Odds:** Multiply the sum of ordinates at odd positions ($y_1, y_3...$) by 4.
5.  **Sum Evens:** Multiply the sum of ordinates at even positions ($y_2, y_4...$) by 2.
6.  **Calculate Total:** Sum all components and multiply by $h/3$.

**Pseudocode**
```text
Input: Function f(x), lower_limit a, upper_limit b, intervals n
Output: Integral value

If n % 2 != 0:
    Print "Error: n must be even"
    Return

h = (b - a) / n
sum = f(a) + f(b) // First and last terms

For i from 1 to n-1:
    x = a + i * h
    If i % 2 == 0:
        sum = sum + 2 * f(x) // Even index
    Else:
        sum = sum + 4 * f(x) // Odd index

Result = sum * (h / 3)
Return Result
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x)
{
    return 1.0 / (1.0 + x * x);
}

double exact_integral(double a, double b)
{
    return atan(b) - atan(a);
}

int main()
{
    ifstream infile("inputsimpson.txt");
    ofstream outfile("outputsimpson13.txt");

    if (!infile || !outfile)
        return 1;

    double a, b;
    int n;
    int caseNum = 1;

    while (infile >> a >> b >> n)
    {
        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Lower Limit: " << a << ", Upper Limit: " << b << ", Subintervals: " << n << endl;

        if (n % 2 != 0)
        {
            outfile << "Error: n must be even for Simpson's 1/3 Rule." << endl;
            outfile << endl;
            continue;
        }

        double h = (b - a) / n;
        double sum = f(a) + f(b);

        for (int i = 1; i < n; i++)
        {
            if (i % 2 == 0)
            {
                sum += 2 * f(a + i * h);
            }
            else
            {
                sum += 4 * f(a + i * h);
            }
        }

        double result = (h / 3.0) * sum;
        double exact = exact_integral(a, b);
        double error = abs(exact - result);

        outfile << "Calculated Value: " << fixed << setprecision(6) << result << endl;
        outfile << "Exact Value:      " << fixed << setprecision(6) << exact << endl;
        outfile << "Absolute Error:   " << scientific << setprecision(4) << error << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
0 1 6
0 6 12
1 2 18
```

#### Output
```
Case 1:
Lower Limit: 0, Upper Limit: 1, Subintervals: 6
Calculated Value: 0.785398
Exact Value:      0.785398
Absolute Error:   2.1816e-007

Case 2:
Lower Limit: 0.0000e+000, Upper Limit: 6.0000e+000, Subintervals: 12
Calculated Value: 1.403702
Exact Value:      1.405648
Absolute Error:   1.9455e-003

Case 3:
Lower Limit: 1.0000e+000, Upper Limit: 2.0000e+000, Subintervals: 18
Calculated Value: 0.321751
Exact Value:      0.321751
Absolute Error:   1.1895e-008


```

---
**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/F.%20Numerical%20Integration/Simpson%201%20by%203)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Simpson’s 1/3 Rule - GeeksforGeeks (Implementation)](https://www.geeksforgeeks.org/program-simpsons-13-rule/)
* [Simpson's Rule Derivation - Wolfram MathWorld](https://mathworld.wolfram.com/SimpsonsRule.html)

<a id="simpson-3-by-8"></a>
### 2. Simpson’s 3/8 Rule

**Theory: Cubic Approximation**
While the 1/3 rule uses parabolas (3 points), Simpson’s 3/8 Rule fits a **third-order polynomial (cubic curve)** through every four points. This generally provides slightly better accuracy for functions that are smoother.



The constraint for this method is that the number of intervals $n$ must be a **multiple of 3**.

$$I \approx \frac{3h}{8} \left[ (y_0 + y_n) + 3(y_1 + y_2 + y_4 + y_5 + \dots) + 2(y_3 + y_6 + \dots) \right]$$

**Algorithm**
1.  **Verify Intervals:** Ensure $n$ is a multiple of 3.
2.  **Calculate Step Size:** $h = (b - a) / n$.
3.  **Sum Extremes:** Add $y_0$ and $y_n$.
4.  **Sum Multiples of 3:** Multiply terms at indices divisible by 3 ($y_3, y_6...$) by 2.
5.  **Sum Others:** Multiply all remaining terms ($y_1, y_2, y_4, y_5...$) by 3.
6.  **Calculate Total:** Sum components and multiply by $3h/8$.

**Pseudocode**
```text
Input: Function f(x), lower_limit a, upper_limit b, intervals n
Output: Integral value

If n % 3 != 0:
    Print "Error: n must be divisible by 3"
    Return

h = (b - a) / n
sum = f(a) + f(b)

For i from 1 to n-1:
    x = a + i * h
    If i % 3 == 0:
        sum = sum + 2 * f(x) // Multiple of 3
    Else:
        sum = sum + 3 * f(x) // Rest of the terms

Result = sum * (3 * h / 8)
Return Result
```
#### Code
```cpp
#include <bits/stdc++.h>

using namespace std;

double f(double x)
{
    return 1.0 / (1.0 + x * x);
}

double exact_integral(double a, double b)
{
    return atan(b) - atan(a);
}

int main()
{
    ifstream infile("inputsimpson.txt");
    ofstream outfile("outputsimpson38.txt");

    if (!infile || !outfile)
        return 1;

    double a, b;
    int n;
    int caseNum = 1;

    while (infile >> a >> b >> n)
    {
        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Lower Limit: " << a << ", Upper Limit: " << b << ", Subintervals: " << n << endl;

        if (n % 3 != 0)
        {
            outfile << "Error: n must be a multiple of 3 for Simpson's 3/8 Rule." << endl;
            outfile << endl;
            continue;
        }

        double h = (b - a) / n;
        double sum = f(a) + f(b);

        for (int i = 1; i < n; i++)
        {
            if (i % 3 == 0)
            {
                sum += 2 * f(a + i * h);
            }
            else
            {
                sum += 3 * f(a + i * h);
            }
        }

        double result = (3.0 * h / 8.0) * sum;
        double exact = exact_integral(a, b);
        double error = abs(exact - result);

        outfile << "Calculated Value: " << fixed << setprecision(3) << result << endl;
        outfile << "Exact Value:      " << fixed << setprecision(3) << exact << endl;
        outfile << "Absolute Error:   " << scientific << setprecision(3) << error << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
0 1 6
0 6 12
1 2 18
```

#### Output
```
Case 1:
Lower Limit: 0, Upper Limit: 1, Subintervals: 6
Calculated Value: 0.785
Exact Value:      0.785
Absolute Error:   2.301e-006

Case 2:
Lower Limit: 0.000e+000, Upper Limit: 6.000e+000, Subintervals: 12
Calculated Value: 1.400
Exact Value:      1.406
Absolute Error:   6.037e-003

Case 3:
Lower Limit: 1.000e+000, Upper Limit: 2.000e+000, Subintervals: 18
Calculated Value: 0.322
Exact Value:      0.322
Absolute Error:   2.609e-008


```

---
**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/F.%20Numerical%20Integration/Simpson%203%20by%208)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Simpson’s 3/8 Rule - GeeksforGeeks](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule)
* [Numerical Integration Rules - Paul's Online Math Notes](https://tutorial.math.lamar.edu/Classes/CalcII/ApproximatingDefIntegrals.aspx)

---

## G. Curve Fitting Methods

<a id="linear"></a>
### 1. Linear Equation 

**Theory: Fitting a Straight Line**
This is the simplest form of regression. We assume the relationship between the dependent variable $y$ and independent variable $x$ is a straight line:
$$y = a_0 + a_1x$$

To find the best $a_0$ (intercept) and $a_1$ (slope), we minimize the error $S = \sum (y_i - (a_0 + a_1x_i))^2$. Taking the partial derivatives with respect to $a_0$ and $a_1$ and setting them to zero gives us the **Normal Equations**:

1.  $\sum y = n \cdot a_0 + a_1 \sum x$
2.  $\sum xy = a_0 \sum x + a_1 \sum x^2$

Solving this system gives direct formulas for the coefficients:
$$a_1 = \frac{n\sum xy - \sum x \sum y}{n\sum x^2 - (\sum x)^2}$$
$$a_0 = \bar{y} - a_1\bar{x}$$
*(Where $\bar{y}$ and $\bar{x}$ are the mean values of y and x, and $n$ is the number of points)*

**Algorithm**
1.  **Initialize Sums:** Set variables for $\sum x$, $\sum y$, $\sum xy$, and $\sum x^2$ to zero.
2.  **Accumulate Data:** Loop through all $n$ data points. For each point $(x_i, y_i)$, add values to the respective sums.
3.  **Calculate Slope ($a_1$):** Apply the formula using the calculated sums.
4.  **Calculate Intercept ($a_0$):** Use the means $\bar{x}$ and $\bar{y}$ and the slope $a_1$.
5.  **Construct Model:** The final equation is $y = a_0 + a_1x$.

**Pseudocode**
```text
Input: Arrays x[] and y[], integer n
Output: Slope a1, Intercept a0

sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0

For i from 0 to n-1:
    sum_x  = sum_x + x[i]
    sum_y  = sum_y + y[i]
    sum_xy = sum_xy + (x[i] * y[i])
    sum_x2 = sum_x2 + (x[i] * x[i])

// Denominator for the slope formula
denom = (n * sum_x2) - (sum_x * sum_x)

If denom == 0:
    Print "Error: Denominator is zero"
    Return

a1 = ((n * sum_xy) - (sum_x * sum_y)) / denom
a0 = (sum_y / n) - (a1 * (sum_x / n))

Print "Equation: y = " + a0 + " + " + a1 + "x"
```
#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream infile("inputlinear.txt");
    ofstream outfile("outputlinear.txt");

    if (!infile || !outfile)
        return 1;

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n), y(n);
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;

        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i];
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            sum_x2 += x[i] * x[i];
        }

        double b = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
        double a = (sum_y - b * sum_x) / n;

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Data Points: " << n << endl;
        outfile << "Linear Fit Equation: y = " << fixed << setprecision(3) << a << " + " << b << "x" << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
4
1 3
2 5
3 7
4 9
5
1 2
2 4.1
3 5.9
4 8.1
5 10
```

#### Output
```
Case 1:
Data Points: 4
Linear Fit Equation: y = 1.000 + 2.000x

Case 2:
Data Points: 5
Linear Fit Equation: y = 0.020 + 2.000x


```

---
**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/G.%20Curve%20Fitting/Linear)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Linear Regression - Yale University](http://www.stat.yale.edu/Courses/1997-98/101/linreg.htm)
* [Least Squares Regression - MathWorld](https://mathworld.wolfram.com/LeastSquaresFitting.html)

<a id="polynomial"></a>
### 2. Polynomial Equation

**Theory: Extending to Higher Orders**
When data shows a curve with peaks and valleys, a straight line is insufficient. We can fit a polynomial of degree $m$:
$$y = a_0 + a_1x + a_2x^2 + \dots + a_mx^m$$

Although this looks non-linear in terms of $x$, it is **linear in terms of the coefficients** $a_0, a_1, \dots$. We minimize the squared error just like before. This results in a system of $(m+1)$ linear equations (Normal Equations) that must be solved simultaneously (often using Gaussian Elimination).

For a 2nd-degree polynomial (Parabola: $y = a_0 + a_1x + a_2x^2$), the matrix system is:

```math
\begin{bmatrix}
n & \sum x & \sum x^2 \\
\sum x & \sum x^2 & \sum x^3 \\
\sum x^2 & \sum x^3 & \sum x^4
\end{bmatrix}
\begin{bmatrix}
a_0 \\
a_1 \\
a_2
\end{bmatrix}
=
\begin{bmatrix}
\sum y \\
\sum xy \\
\sum x^2y
\end{bmatrix}
```
**Algorithm**
1.  **Select Degree $m$:** Choose the order of the polynomial (e.g., $m=2$ for a parabola).
2.  **Calculate Power Sums:** Compute sums for $x$ up to power $2m$ ($\sum x, \sum x^2, \dots \sum x^{2m}$).
3.  **Calculate Moment Sums:** Compute sums like $\sum y, \sum xy, \sum x^2y$.
4.  **Build Matrix:** Fill the augmented matrix with these sums.
5.  **Solve System:** Use Gaussian Elimination to solve for $a_0, a_1, \dots a_m$.

**Pseudocode (Building the Matrix for Degree m)**
```text
Input: Arrays x[] and y[], degree m, points n
Output: Coefficients a[]

// 1. Calculate the sums of powers of x
powers[] = array of size (2*m + 1)
For k from 0 to 2*m:
    powers[k] = sum of (x[i]^k) for all i

// 2. Calculate sums of y * x^k
moments[] = array of size (m + 1)
For k from 0 to m:
    moments[k] = sum of (y[i] * x[i]^k) for all i

// 3. Construct Matrix System (Matrix B, Vector C)
For i from 0 to m:
    For j from 0 to m:
        B[i][j] = powers[i + j]
    C[i] = moments[i]

// 4. Solve B * a = C
coefficients = GaussianElimination(B, C)
Return coefficients
```
#### Code
```cpp
#include <bits/stdc++.h>

using namespace std;

void gaussianElimination(vector<vector<double>> &a, vector<double> &b, vector<double> &x, int n)
{
    for (int i = 0; i < n; i++)
    {
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (abs(a[k][i]) > abs(a[maxRow][i]))
            {
                maxRow = k;
            }
        }

        swap(a[i], a[maxRow]);
        swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; k++)
        {
            double factor = a[k][i] / a[i][i];
            for (int j = i; j < n; j++)
            {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}

int main()
{
    ifstream infile("inputpoly.txt");
    ofstream outfile("outputpoly.txt");

    if (!infile || !outfile)
        return 1;

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i];
        }

        int degree;
        infile >> degree;
        int m = degree + 1;

        vector<vector<double>> A(m, vector<double>(m, 0.0));
        vector<double> B(m, 0.0);
        vector<double> coeffs(m);

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                {
                    sum += pow(x[k], i + j);
                }
                A[i][j] = sum;
            }
            double sum = 0;
            for (int k = 0; k < n; k++)
            {
                sum += y[k] * pow(x[k], i);
            }
            B[i] = sum;
        }

        gaussianElimination(A, B, coeffs, m);

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Data Points: " << n << ", Degree: " << degree << endl;
        outfile << "Polynomial Fit Equation: y = ";
        for (int i = 0; i < m; i++)
        {
            if (i == 0)
                outfile << fixed << setprecision(4) << coeffs[i];
            else if (coeffs[i] >= 0)
                outfile << " + " << coeffs[i] << "x^" << i;
            else
                outfile << " - " << abs(coeffs[i]) << "x^" << i;
        }
        outfile << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
4
0 1
1 1.8
2 1.3
3 2.5
2
5
0 1
1 6
2 17
3 34
4 57
2
```

#### Output
```
Case 1:
Data Points: 4, Degree: 2
Polynomial Fit Equation: y = 1.1500 + 0.1000x^1 + 0.1000x^2

Case 2:
Data Points: 5, Degree: 2
Polynomial Fit Equation: y = 1.0000 + 2.0000x^1 + 3.0000x^2


```

---

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/G.%20Curve%20Fitting/Polynomial)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Polynomial Regression - GeeksforGeeks](https://www.geeksforgeeks.org/polynomial-regression-for-non-linear-data-ml/)
* [Least Squares Fitting of Polynomials - Wolfram MathWorld](https://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html)

<a id="transcendental"></a>
### 3. Transcendental Equation

**Theory: Linearizatio of Non-Linear Models**
Sometimes data does not fit a straight line but follows an exponential ($y = ae^{bx}$) or power ($y = ax^b$) law. We cannot apply the standard least-squares formulas directly to these non-linear forms.



Instead, we **linearize** the equation using logarithms:
* **Exponential Model ($y = ae^{bx}$):** Take $\ln$ of both sides $\rightarrow \ln(y) = \ln(a) + bx$.
    * This looks like a line $Y = A_0 + A_1x$, where $Y = \ln(y)$, $A_0 = \ln(a)$, and $A_1 = b$.
* **Power Model ($y = ax^b$):** Take $\log$ of both sides $\rightarrow \log(y) = \log(a) + b\log(x)$.
    * This looks like $Y = A_0 + A_1X$, where $Y = \log(y)$ and $X = \log(x)$.

Once linearized, we calculate the slope and intercept using the standard linear formulas, then transform them back to find the original constants $a$ and $b$.

**Algorithm (for Exponential $y = ae^{bx}$)**
1.  **Transform Data:** Create a new array $Z$ where $Z_i = \ln(y_i)$.
2.  **Apply Linear Regression:** Perform standard linear regression on the pairs $(x_i, Z_i)$.
3.  **Calculate $a_1$ (slope):** This corresponds directly to $b$.
4.  **Calculate $a_0$ (intercept):** This corresponds to $\ln(a)$.
5.  **Inverse Transform:** Calculate $a = e^{a_0}$.
6.  **Final Model:** $y = a e^{bx}$.

**Pseudocode**
```text
Input: Arrays x[] and y[], integer n
Output: Coefficients a and b for y = ae^(bx)

sum_x = 0, sum_z = 0, sum_xz = 0, sum_x2 = 0

For i from 0 to n-1:
    z = ln(y[i])  // Linearize y
    sum_x  = sum_x + x[i]
    sum_z  = sum_z + z
    sum_xz = sum_xz + (x[i] * z)
    sum_x2 = sum_x2 + (x[i] * x[i])

denom = (n * sum_x2) - (sum_x * sum_x)
b = ((n * sum_xz) - (sum_x * sum_z)) / denom  // Slope
A0 = (sum_z / n) - (b * (sum_x / n))          // Intercept

a = exp(A0) // Inverse transform

Print "Equation: y = " + a + " * e^(" + b + "x)"
```
#### Code
```cpp
#include <bits/stdc++.h>

using namespace std;

int main()
{
    ifstream infile("inputtrans.txt");
    ofstream outfile("outputtrans.txt");

    if (!infile || !outfile)
        return 1;

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n), y(n);
        double sum_x = 0, sum_Y = 0, sum_xY = 0, sum_x2 = 0;

        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i];
            double Y = log(y[i]);
            sum_x += x[i];
            sum_Y += Y;
            sum_xY += x[i] * Y;
            sum_x2 += x[i] * x[i];
        }

        double B = (n * sum_xY - sum_x * sum_Y) / (n * sum_x2 - sum_x * sum_x);
        double A = (sum_Y - B * sum_x) / n;

        double a = exp(A);
        double b = B;

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Data Points: " << n << endl;
        outfile << "Transcendental Fit Equation (y = ae^bx):" << endl;
        outfile << "y = " << fixed << setprecision(4) << a << "e^(" << b << "x)" << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}

```

#### Input
```
5
1 0.5
2 2.0
3 4.5
4 8.0
5 12.5
4
1 2.7
2 7.4
3 20.1
4 54.6
```

#### Output
```
Case 1:
Data Points: 5
Transcendental Fit Equation (y = ae^bx):
y = 0.3245e^(0.7824x)

Case 2:
Data Points: 4
Transcendental Fit Equation (y = ae^bx):
y = 0.9940e^(1.0020x)


```

---
**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/G.%20Curve%20Fitting/Transcendental)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Curve Fitting - Wikipedia](https://en.wikipedia.org/wiki/Curve_fitting)
* [Exponential Curve Fitting - MathWorld](https://mathworld.wolfram.com/LeastSquaresFittingExponential.html)

---

## Sample Run

### Input

```text
# Give Input in the input file 
Enter number of equations: 3
Enter coefficients...
```

### Output

```text
Solution:
x = 1.0
y = 2.0
z = 3.0
```

---

## Contribution

### Team Details

| Full Name | GitHub Username | Roll Number |
| :--- | :--- | :--- |
| **MD SUAIB AHMED SAFI** | [@suaib022](https://github.com/suaib022) | 2207115 |
| **ASHRAFUR RAHMAN NIHAD** | [@ARN101](https://github.com/ARN101) | 2207116 |
| **DADHICHI SAREKR SHAYON** | [@Dadhichi-Sarker-Shayon](https://github.com/Dadhichi-Sarker-Shayon) | 2207118 |

---

### Individual Technical Contributions

**MD SUAIB AHMED SAFI:**
Responsible for implementing the LU Decomposition Method, Runge–Kutta Method, and Matrix Inverse. Additionally, he handled the logic for Newton Forward Interpolation, Newton Backward Interpolation, and Linear Equation solvers.

**ASHRAFUR RAHMAN NIHAD:**
Developed the code for the Jacobi & Gauss–Seidel Method, Simpson’s 1/3 Rule, and Simpson’s 3/8 Rule. He also contributed the implementations for Newton Forward Differentiation, Polynomial Equation, and Transcendental Equation.

**DADHICHI SAREKR SHAYON:**
Focused on the direct linear algebra methods including the Gauss Elimination Method and Gauss Jordan Method. His contributions also cover root-finding algorithms such as the Bisection Method, False Position Method, Secant Method, and Newton Raphson Method.


---

## 🤝 Future Contribution Guidelines

We welcome contributions to improve this project! Whether you want to fix a bug, add a new numerical method, or improve documentation, here is how you can help.

### How to Contribute

1.  **Fork the Repository**
    - Click the "Fork" button at the top right of this page to create your own copy of the repository.

2.  **Clone Your Fork**
    ```bash
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    cd your-repo-name
    ```

3.  **Create a New Branch**
    - Always create a separate branch for your changes. Do not push directly to `main`.
    ```bash
    git checkout -b feature/your-feature-name
    # OR for bug fixes
    git checkout -b fix/bug-name
    ```

4.  **Make Your Changes**
    - Ensure your code follows the existing style (indentation, variable naming, comments).
    - If adding a new method, please include a sample input/output test case in the comments.

5.  **Commit Your Changes**
    - Write clear, concise commit messages.
    ```bash
    git add .
    git commit -m "Added implementation for [Method Name]"
    ```

6.  **Push to Your Fork**
    ```bash
    git push origin feature/your-feature-name
    ```

7.  **Submit a Pull Request (PR)**
    - Go to the original repository and click **"Compare & pull request"**.
    - Provide a brief description of what you added or fixed.
    - Tag the maintainers for review.

### Code Standards
- **Language:** C++ (Standard 11 or higher) / Python (if applicable).
- **Documentation:** Every function must have comments explaining the parameters and return values.
- **Testing:** Verify that your algorithm works on edge cases (e.g., division by zero checks).

---

### 🐛 Reporting Issues
Found a bug? Please open an **Issue** with the following details:
- The method causing the error.
- Input data used.
- Expected vs. Actual output.
---

## References

1. Chapra, S. C., & Canale, R. P. *Numerical Methods for Engineers*
2. Burden & Faires, *Numerical Analysis*
3. NPTEL Lecture Notes on Numerical Methods
4. *Numerical Methods for Engineers* — Steven Chapra
5. MIT OpenCourseWare (Numerical Analysis)
6. NPTEL Numerical Methods


