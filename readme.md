# Numerical Methods Console Application

## Overview

This project is a comprehensive collection of numerical methods implemented in C++. It covers fundamental algorithms for solving linear and non-linear equations, interpolation, numerical differentiation and integration, differential equations, and curve fitting techniques.

---

## Project Structure

```text
nm/
├── readme.md
│
├── A. Solution of Linear Equations/
│   ├── Gauss Elimination/
│   ├── Gauss Jordan/
│   ├── LU Decomposition/
│   ├── Matrix Inverse/
│   └── Iterative Methods/
│       ├── Jacobi/
│       └── GaussSeidel/
│
├── B. Solution of Non-Linear Equations/
│   ├── Bisection/
│   ├── False_Position/
│   ├── Secant/
│   └── Newton_Raphson/
│
├── C. Interpolation and Approximation/
│   ├── Newton Forward/
│   └── Newton Backward/
│
├── E. Solution of Differential Equations/
│   ├── Newton Forward Differentiation/
│   └── Runge Kutta/
│
├── F. Numerical Integration/
│   ├── Simpson 1 by 3/
│   └── Simpson 3 by 8/
│
└── G. Curve Fitting/
    ├── Linear/
    ├── Polynomial/
    └── Transcendental/
```

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
    * *Note: The system must be diagonally dominant ($|a_{ii}| > \sum |a_{ij}|$) for guaranteed convergence.*
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
The Bisection Method is a simple and dynamic numerical technique used to find roots of continuous functions. It is also known as:
1. Binary chopping method
2. Half-interval method

**Algorithm of Bisection Method:**

Step 1: Choosing two numbers (x1) and (x2) such that:
[ f(x1) f(x2) < 0 ]

Step 2: Computing the midpoint:
[ x0 = (x1 + x2)/2 ]

Step 3: Evaluating f(x0).

Step 4:
If (f(x0) == 0), then (x0) is the exact root.
If (f(x0) .f(x1) < 0), setting (x2 = x0).
If (f(x0) .f(x2) < 0), setting (x1 = x0).

Step 5: Repeating until the stopping criterion is met:
|x2 - x1| < E

Step 6: Stopping. The approximate root is (x0).

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

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/B.%20Solution%20of%20Non-Linear%20Equations/Bisection)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Bisection Method - Wikipedia](https://en.wikipedia.org/wiki/Bisection_method)
- [Bisection Method - GeeksforGeeks](https://www.geeksforgeeks.org/bisection-method/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="false_position"></a>
### 2. False Position Method

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
The Secant Method approximates the derivative numerically using two points instead of requiring an analytical derivative. Given two initial points x(n-1) and xn :

f'(xn) = (f(xn) - f(x(n-1))) / (xn - x(n-1))

Substituting this into Newton's formula:

[x(n+1) = xn - f(xn) * (xn - x(n-1)) / (f(xn) - f(x(n-1)))]

**Algorithm of Secant Method**

Step 1: Choosing two initial approximations x0 and x1.
Step 2: Applying:
[ x(n+1) = xn - f(xn) * (xn - x(n-1)) / (f(xn) - f(x(n-1))) ]
Step 3: Checking:
| x(n+1) - xn | < E

If yes -> Stop.
Else -> set x(n-1) = xn, xn = x(n+1) and repeat.

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
The Newton-Raphson method uses the tangent line at the current approximation to estimate a better root approximation. Given a guess ( xn ), the tangent line at that point is:

[y = f(xn) + f'(xn)(x - xn)]

Setting (y = 0): [0 = f(xn) + f'(xn)(x - xn)]

Solving for (x): [x(n+1) = xn - f(xn)/f'(xn)]

This is the Newton-Raphson formula.

**Algorithm of Newton-Raphson Method**

Step 1: Choosing an initial guess x0.
Step 2: Computing derivative f'(xn).
Step 3: Computing next approximation:
    [x(n+1) = xn - f(xn)/f'(xn)]
Step 4: Checking stopping criterion:
    |x(n+1)-xn| < E

If satisfied => Stop.
Else => setting (xn = x(n+1)) and repeating.

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

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/C.%20Interpolation%20and%20Approximation/Newton%20Forward)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Newton Forward Interpolation - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- [Newton Forward Interpolation - GeeksforGeeks](https://www.geeksforgeeks.org/newton-forward-difference-interpolation/)
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

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/C.%20Interpolation%20and%20Approximation/Newton%20Backward)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
- [Newton Backward Interpolation - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- [Newton Backward Interpolation - GeeksforGeeks](https://www.geeksforgeeks.org/newton-backward-difference-interpolation/)
- Numerical Methods for Engineers - Chapra & Canale

---


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

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/F.%20Numerical%20Integration/Simpson%203%20by%208)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Simpson’s 3/8 Rule - GeeksforGeeks](https://www.geeksforgeeks.org/simpsons-38-rule-python/)
* [Numerical Integration Rules - Swarthmore College](https://lpsa.swarthmore.edu/NumInt/NumIntMain.html)

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

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/G.%20Curve%20Fitting/Transcendental)

Includes:
- C++ source code
- Input file
- Output file

**Further Study**
* [Curve Fitting: Exponential and Power Laws - Math For Engineers](https://www.google.com/search?q=https://www.mathforcollege.com/nm/topics/textbook_index.html)
* [Linearization of Exponential Models - Ximera OSU](https://www.google.com/search?q=https://ximera.osu.edu/mooculus/calculus2/linearization/digInLinearization)

---

## Environment Setup

### Requirements
- C++ Compiler (GCC / MinGW / Clang)
- C++11 or later

### Compile & Run

```bash
g++ code.cpp -o run
./run
```

---

## Sample Run

### Input

```text
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

### Our Team Contribution

| Full Name | GitHub Username | Roll Number |
|-----------|-----------------|-------------|
| MD SUAIB AHMED SAFI | [suaib022](https://github.com/suaib022) | 2207115 |
| ASHRAFUR RAHMAN NIHAD | [ARN101](https://github.com/ARN101) | 2207116 |
| DADHICHI SAREKR SHAYON | [Dadhichi-Sarker-Shayon](https://github.com/Dadhichi-Sarker-Shayon) | 2207118 |

- **Member 1:** Method implementation & testing
- **Member 2:** Mathematical theory & PDFs
- **Member 3:** Repository structure & documentation

---

## Future Contributions

We welcome contributions! 🙌

**How to contribute:**
1. Fork the repository
2. Create a new branch
3. Add method / optimize code / improve documentation
4. Submit a pull request

---

## References

1. Chapra, S. C., & Canale, R. P. *Numerical Methods for Engineers*
2. Burden & Faires, *Numerical Analysis*
3. NPTEL Lecture Notes on Numerical Methods
4. *Numerical Methods for Engineers* — Steven Chapra
5. MIT OpenCourseWare (Numerical Analysis)
6. NPTEL Numerical Methods
