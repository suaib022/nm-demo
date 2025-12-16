# Numerical Methods

## Summary
### A. Solution of Linear Equations
*   Gauss Elimination
*   Gauss-Jordan Elimination
*   [LU Factorization](#lu-decomposition-method)
*   Iterative Methods: Jacobi and Gauss-Seidel methods

### B. Solution of Non-linear Equations
*   Bisection Method
*   False Position Method (Regula-Falsi)
*   Secant Method
*   Newton-Raphson Method

### C. Interpolation and Approximation
*   [Newton Forward Interpolation](#newton-forward-interpolation)
*   [Newton Backward Interpolation](#newton-backward-interpolation)
*   Error Analysis
*   [Newton Divided Difference Interpolation](#newton-divided-difference-interpolation)

### D. Numerical Differentiation
*   Equal-Interval Interpolation Method
*   Second-Order Derivative Formula
*   Lagrange’s Interpolation-Based Differentiation

### E. Solution of Differential Equations
*   [Runge-Kutta Method](#runge-kutta-rk-method)

### F. Numerical Integration
*   Simpson’s 1/3 Rule
*   Simpson’s 3/8 Rule
*   Trapezium/Trapezoidal Rule
*   One-third Rule (alternate naming for Simpson’s 1/3)

### G. Curve Fitting
*   Least-Squares Straight Lines
*   Least-Squares Polynomials
*   Non-Linear Curve Fitting

### H. Matrix Inversion Approach

---

## A. Solution of Linear Equations

In numerical analysis, solving linear systems (typically written as $Ax = b$) is fundamental. While you might be familiar with "Direct Methods" like Cramer's Rule or Gaussian Elimination that attempt to find the exact solution in a finite number of steps, they can become computationally expensive for very large systems.

This is where **Iterative Methods** shine. Instead of trying to solve the problem in one go, these methods start with a guess and refine it over and over again until the error is negligible.

### 1. Iterative Methods: The Art of Successive Refinement

**Why "Iterative"?**
The term comes from the Latin *iterare* (to repeat). Unlike direct elimination, these algorithms generate a sequence of approximate solutions $\{x^{(0)}, x^{(1)}, x^{(2)}, ...\}$. Each step "iterates" on the previous one to reduce the error. Think of it like tuning a guitar: you pluck the string, check the pitch, adjust the peg slightly, and repeat until it sounds perfect.

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
**Further Study**
* [Jacobi Method - Wikipedia (Comprehensive Theory)](https://en.wikipedia.org/wiki/Jacobi_method)
* [Jacobi Method Explained - GeeksforGeeks (Examples & Code)](https://www.geeksforgeeks.org/engineering-mathematics/jacobian-method/)

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
**Further Study**
* [Gauss-Seidel Method - GeeksforGeeks](https://www.geeksforgeeks.org/gauss-seidel-method/)
* [Iterative Methods: Jacobi vs Gauss-Seidel - Math LibreTexts](https://math.libretexts.org/Bookshelves/Linear_Algebra/Introduction_to_Matrix_Algebra_(Kaw)/01%3A_Chapters/1.08%3A_Gauss-Seidel_Method)

<a id="lu-decomposition-method"></a>
### LU Decomposition Method

**Introduction**
LU decomposition is used in numerical methods to solve linear equations. To solve a problem like **Ax = B**, LU decomposition breaks the matrix (**A**) into two simpler matrices: a lower-triangular matrix (**L**) and an upper-triangular matrix (**U**).

Another big advantage is efficiency: if we need to solve the same matrix (**A**) with many different **b** vectors, we only need to do LU factorization once.

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

---

## C. Interpolation and Approximation

<a id="newton-forward-interpolation"></a>
### Newton Forward Interpolation

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

<a id="newton-backward-interpolation"></a>
### Newton Backward Interpolation

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

<a id="newton-divided-difference-interpolation"></a>
### Newton Divided Difference Interpolation

**Theory**
Newton's Divided Difference Interpolation is a method for constructing an interpolating polynomial for a given set of data points where the interval between data points is **not necessarily equal**.

The formula is given by:
$$f(x) = f(x_0) + (x-x_0)f[x_0, x_1] + (x-x_0)(x-x_1)f[x_0, x_1, x_2] + \dots + (x-x_0)\dots(x-x_{n-1})f[x_0, \dots, x_n]$$

**Algorithm**
1.  **Input**: Read ($n$) data points ($x, y$) and the value to interpolate ($value$).
2.  **Table**: Construct the divided difference table using the recursive formula.
3.  **Compute**: Apply formula.
4.  **Output**: Display the interpolated value.

**Pseudocode**
```cpp
function newton_divided(x[], y[][], n, value):
    for i from 1 to n-1:
        for j from 0 to n-i-1:
            y[j][i] = (y[j+1][i-1] - y[j][i-1]) / (x[i+j] - x[j])

    sum = y[0][0]
    for i from 1 to n-1:
        pro = 1
        for j from 0 to i-1:
            pro = pro * (value - x[j])
        sum = sum + (pro * y[0][i])

    return sum
```

---

## D. Numerical Differentiation

Numerical differentiation is the process of calculating the derivative (rate of change) of a function using a set of discrete data points $(x_i, y_i)$ rather than an analytical formula.

### 1. Equal-Interval Interpolation Method

**Theory: Differentiating the Polynomial**
When data points are spaced equally, we can use Newton's Forward Difference formula.
$$\frac{dy}{dx} = \frac{1}{h} \left[ \Delta y_0 + \frac{2p-1}{2} \Delta^2 y_0 + \frac{3p^2-6p+2}{6} \Delta^3 y_0 + \dots \right]$$

**Algorithm**
1.  **Check Interval:** Verify that all $x$ values have a constant difference $h$.
2.  **Difference Table:** Construct the Forward Difference Table.
3.  **Calculate $p$:** Determine the position factor.
4.  **Apply Series:** Substitute the difference values.
5.  **Scale:** Divide the result by step size $h$.

**Pseudocode**
```text
Input: Arrays x[] and y[], value xp (point to differentiate at)
Output: Derivative value dy/dx

1. Calculate step size h = x[1] - x[0]
2. Generate Forward Difference Table diff[n][n]
   For j = 1 to n-1:
       For i = 0 to n-j-1:
           diff[i][j] = diff[i+1][j-1] - diff[i][j-1]

3. Find index i such that x[i] is closest to xp
4. Calculate p = (xp - x[i]) / h
5. Initialize sum = 0

6. Apply Formula (Example for 1st derivative):
   term1 = diff[i][0]
   term2 = (2*p - 1) * diff[i][1] / 2
   term3 = (3*p*p - 6*p + 2) * diff[i][2] / 6
   sum = term1 + term2 + term3 + ...

7. Result = sum / h
8. Return Result
```

---

## E. Solution of Differential Equations

<a id="runge-kutta-rk-method"></a>
### Runge-Kutta (RK) Method

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

---

## F. Numerical Integration

Numerical integration calculates the approximate value of a definite integral $\int_{a}^{b} f(x) dx$.

### 1. Simpson’s 1/3 Rule

**Theory: Parabolic Approximation**
Approximates the function $f(x)$ as a second-order polynomial (parabola). Requires even intervals.
$$I \approx \frac{h}{3} \left[ (y_0 + y_n) + 4(y_1 + y_3 + \dots) + 2(y_2 + y_4 + \dots) \right]$$

**Algorithm**
1.  **Verify Intervals:** Ensure $n$ is even.
2.  **Calculate Step Size:** $h = (b - a) / n$.
3.  **Sum**: Apply formula.

**Pseudocode**
```text
Input: Function f(x), lower_limit a, upper_limit b, intervals n
Output: Integral value

If n % 2 != 0: Returns Error

h = (b - a) / n
sum = f(a) + f(b)

For i from 1 to n-1:
    x = a + i * h
    If i % 2 == 0: sum = sum + 2 * f(x)
    Else: sum = sum + 4 * f(x)

Result = sum * (h / 3)
Return Result
```

### 2. Simpson’s 3/8 Rule

**Theory: Cubic Approximation**
Fits a third-order polynomial. Requires intervals to be a multiple of 3.
$$I \approx \frac{3h}{8} \left[ (y_0 + y_n) + 3(y_1 + y_2 + y_4 + \dots) + 2(y_3 + y_6 + \dots) \right]$$

**Algorithm**
1.  **Verify Intervals:** Ensure $n$ % 3 == 0.
2.  **Calculate Step Size**: $h$.
3.  **Sum**: Apply formula.

**Pseudocode**
```text
Input: Function f(x), lower_limit a, upper_limit b, intervals n
Output: Integral value

If n % 3 != 0: Return Error

h = (b - a) / n
sum = f(a) + f(b)

For i from 1 to n-1:
    x = a + i * h
    If i % 3 == 0: sum = sum + 2 * f(x)
    Else: sum = sum + 3 * f(x)

Result = sum * (3 * h / 8)
Return Result
```

---

## G. Curve Fitting (Regression)

Curve fitting constructs a mathematical function that best fits a series of data points, minimizing error (Least Squares).

### 1. Least-Squares Regression: Linear Equation
**Theory**: Fits $y = a_0 + a_1x$ by solving Normal Equations.
**Pseudocode**:
```text
Input: Arrays x[] and y[], integer n
sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0
// ... Calculate sums ...
denom = (n * sum_x2) - (sum_x * sum_x)
a1 = ((n * sum_xy) - (sum_x * sum_y)) / denom
a0 = (sum_y / n) - (a1 * (sum_x / n))
Print "y = " + a0 + " + " + a1 + "x"
```

### 2. Least-Squares Regression: Transcendental Equation
**Theory**: Linearize $y = ae^{bx}$ to $\ln(y) = \ln(a) + bx$.
**Pseudocode**:
```text
Input: Arrays x[] and y[], integer n
// ... Calculate sums with z = ln(y) ...
b = ((n * sum_xz) - (sum_x * sum_z)) / denom
A0 = (sum_z / n) - (b * (sum_x / n))
a = exp(A0)
Print "y = " + a + " * e^(" + b + "x)"
```

### 3. Least-Squares Regression: Polynomial Equation
**Theory**: Fits $y = a_0 + a_1x + \dots + a_mx^m$ by solving a matrix system.
**Pseudocode**:
```text
Input: Arrays x[] and y[], degree m, points n
// 1. Calculate sums of powers of x
// 2. Calculate sums of moments (y * x^k)
// 3. Build Matrix B and Vector C
// 4. Solve B * a = C
Return coefficients
```
