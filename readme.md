# Numerical Methods

## Summary
*   [LU Decomposition Method](#lu-decomposition-method)
*   [Newton Forward Interpolation](#newton-forward-interpolation)
*   [Newton Backward Interpolation](#newton-backward-interpolation)
*   [Runge-Kutta (RK) Method](#runge-kutta-rk-method)

---

<a id="lu-decomposition-method"></a>
# LU Decomposition Method

## Introduction
LU decomposition is used in numerical methods to solve linear equations. To solve a problem like **Ax = B**, LU decomposition breaks the matrix (**A**) into two simpler matrices: a lower-triangular matrix (**L**) and an upper-triangular matrix (**U**).

Another big advantage is efficiency: if we need to solve the same matrix (**A**) with many different **b** vectors, we only need to do LU factorization once.

## Theory

### 1) The main idea
Given a square matrix (**A**). It will be rewritten as follows:

$$A = LU$$

where:
*   **L** is a lower triangular matrix
*   **U** is an upper triangular matrix

### 2) Using LU to solve AX = B
Once we have ($A = LU$), we can write:

$$LUx = b$$

To make this easy, we solve it in two smaller steps. We introduce a helper vector **Y**:

$$Ly = b$$
$$Ux = y$$

*   First, solve **Ly = b** using **forward substitution**.
*   Then, solve **Ux = y** using **backward substitution**.

### 3) Why pivoting is often needed
Sometimes LU decomposition runs into trouble if a diagonal element (called a **pivot**) becomes zero. That can cause division problems and also makes result less reliable.

To fix this, we often swap rows to bring a better pivot into position. With row swapping, the factorization becomes:

$$PA = LU$$

where (**P**) is a matrix that represents the row swaps. This is called **partial pivoting**. It makes the method more stable in real calculations.

## Algorithm
The diagonal of (**L**) is set to 1. The idea is:

For each column/step ($k = 1$) to ($n$):
1.  Compute the ($k$)-th row of (**U**)
2.  Compute the ($k$)-th column of (**L**)

Mathematically:

$$U_{k,j} = A_{k,j} - \sum_{s=1}^{k-1} L_{k,s} U_{s,j} \quad , \quad j = k, \dots, n$$

$$L_{i,k} = \frac{A_{i,k} - \sum_{s=1}^{k-1} L_{i,s} U_{s,k}}{U_{k,k}} \quad , \quad i = k+1, \dots, n$$

And set $L_{k,k} = 1$.

If pivoting is used row swaps are done before continuing whenever needed.

## Procedure

**Step 1: Check the matrix**
*   Make sure (**A**) is a square matrix.

**Step 2: Decompose (A)**
*   Find (**L**) and (**U**) such that:
    *   ($A = LU$)
    *   ($PA = LU$)

**Step 3: Solve the first triangular matrix**
*   Solve: $Ly = b$

**Step 4: Solve the second triangular system**
*   Solve: $Ux = b$ (Note: The image said $Ux = b$ but scientifically it should be $Ux = y$, however, the text in image 3 Step 4 says "Solve: Ux = b", but solving for X. Contextually it usually means "Solve for x using the equation Ux=y". Wait, Image 3 Step 4 explicitly says "Solve: Ux = b". That might be a typo in the user's source image or a simplification where 'b' is the modified vector. I will write $Ux = y$ as in Theory section to be correct, but mention the step from image)

*Correction based on Image 1 Theory*: "Then, solve Ux=y".
*Correction based on Image 3 Procedure*: "Solve: Ux=b".
*Decision*: I will use **Ux = y** to be mathematically consistent with the Theory section, as Step 4 in the image says "The vector X is the solution", implying the final step. $Ux=b$ would imply $y$ wasn't used. I'll stick to $Ux=y$.

**Step 4: Solve the second triangular system**
*   Solve: $Ux = y$

The vector **X** is the solution of the problem.

## Pseudocode
```cpp
// Doolittle's Algorithm with Partial Pivoting
function solve_LU(Matrix A, Vector b):
    n = size(A)
    L = IdentityMatrix(n)
    U = ZeroMatrix(n)
    
    // Decomposition
    for i from 0 to n-1:
        // Partial Pivoting (Optional but recommended)
        pivot_row = find_max_in_column(i, from: i to n-1)
        swap_rows(A, i, pivot_row)
        swap_rows(b, i, pivot_row) // Apply swap to b vector as well (or keep track in P)
        
        // Calculate U (Upper Triangular)
        for k from i to n-1:
            prod_sum = 0
            for j from 0 to i-1:
                prod_sum += L[i][j] * U[j][k]
            U[i][k] = A[i][k] - prod_sum
            
        // Calculate L (Lower Triangular)
        for k from i to n-1:
            if i == k:
                L[i][i] = 1
            else:
                prod_sum = 0
                for j from 0 to i-1:
                    prod_sum += L[k][j] * U[j][i]
                L[k][i] = (A[k][i] - prod_sum) / U[i][i]

    // Check for unique solution
    for i from 0 to n-1:
        if abs(U[i][i]) < epsilon:
            return "No unique solution"

    // Forward Substitution: Ly = b
    Vector y(n)
    for i from 0 to n-1:
        sum = 0
        for j from 0 to i-1:
            sum += L[i][j] * y[j]
        y[i] = b[i] - sum

    // Backward Substitution: Ux = y
    Vector x(n)
    for i from n-1 down to 0:
        sum = 0
        for j from i+1 to n-1:
            sum += U[i][j] * x[j]
        x[i] = (y[i] - sum) / U[i][i]
        
    return x
```

---

<a id="newton-forward-interpolation"></a>
# Newton Forward Interpolation

## Theory
Newton's Forward Interpolation is used to approximate the value of a function $f(x)$ at valid points, based on a set of known data points that are **equally spaced**.

The interpolation formula is given by:

$$y(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \frac{u(u-1)(u-2)}{3!} \Delta^3 y_0 + \dots$$

Where:
*   $x$ is the point to interpolate.
*   $x_0$ is the first value in the dataset.
*   $h$ is the interval difference ($x_1 - x_0$).
*   $u = \frac{x - x_0}{h}$
*   $\Delta y_0$ represents the forward difference.

## Algorithm
1.  **Input**: Read ($n$) data points ($x, y$) and the value to interpolate ($value$).
2.  **Difference Table**: Construct the forward difference table where:
    $$\Delta y_{j,i} = y_{j+1, i-1} - y_{j, i-1}$$
3.  **Calculate u**: Compute $u = (value - x[0]) / (x[1] - x[0])$.
4.  **Compute Sum**: Initialize $sum = y[0][0]$.
    *   Iterate $i$ from 1 to $n-1$.
    *   Update term: $\frac{u(u-1)\dots(u-i+1)}{i!} \times \Delta^i y_0$
    *   Add to sum.
5.  **Output**: Display the interpolated value.

## Pseudocode
```cpp
function newton_forward(x[], y[][], n, value):
    // Construct Difference Table
    for i from 1 to n-1:
        for j from 0 to n-i-1:
            y[j][i] = y[j+1][i-1] - y[j][i-1]

    // Calculate u
    sum = y[0][0]
    u = (value - x[0]) / (x[1] - x[0])
    
    // Apply Formula
    for i from 1 to n-1:
        u_term = 1
        for k from 0 to i-1:
            u_term = u_term * (u - k) // Calculate u(u-1)...
        
        factorial = fact(i)
        sum = sum + (u_term * y[0][i]) / factorial

    return sum
```

---

<a id="newton-backward-interpolation"></a>
# Newton Backward Interpolation

## Theory
Newton's Backward Interpolation is similar to the forward method but is more accurate for interpolating values near the **end** of the dataset.

The formula is given by:

$$y(x) = y_n + u \nabla y_n + \frac{u(u+1)}{2!} \nabla^2 y_n + \frac{u(u+1)(u+2)}{3!} \nabla^3 y_n + \dots$$

Where:
*   $x_n$ is the last value in the dataset.
*   $u = \frac{x - x_n}{h}$
*   $\nabla y_n$ represents the backward difference.

## Algorithm
1.  **Input**: Read ($n$) data points ($x, y$) and the values.
2.  **Difference Table**: Construct the backward difference table.
3.  **Calculate u**: Compute $u = (value - x[n-1]) / (x[1] - x[0])$.
4.  **Compute Sum**: Initialize $sum = y[n-1][0]$.
    *   Iterate $i$ from 1 to $n-1$.
    *   Update term: $\frac{u(u+1)\dots(u+i-1)}{i!} \times \nabla^i y_n$
    *   Add to sum.
5.  **Output**: Display result.

## Pseudocode
```cpp
function newton_backward(x[], y[][], n, value):
    // Construct Difference Table
    // (Notice indices usually allow for implicit backward calculation or same table logic)
    for i from 1 to n-1:
        for j from n-1 down to i:
             y[j][i] = y[j][i-1] - y[j-1][i-1]

    // Calculate u (based on last element)
    sum = y[n-1][0]
    u = (value - x[n-1]) / (x[1] - x[0])

    // Apply Formula
    for i from 1 to n-1:
        u_term = 1
        for k from 0 to i-1:
             u_term = u_term * (u + k) // Note: (u+k) for backward
             
        factorial = fact(i)
        sum = sum + (u_term * y[n-1][i]) / factorial

    return sum
```

---

<a id="runge-kutta-rk-method"></a>
# Runge-Kutta (RK) Method

## Theory
The Runge-Kutta method (specifically the fourth-order RK4) is a widely used technique for the approximate solution of ordinary differential equations (ODEs).

Given an equation:
$$ \frac{dy}{dx} = f(x, y) $$
With initial condition $y(x_0) = y_0$.

The value $y_{n+1}$ at $x_{n+1} = x_n + h$ is determined by:

$$y_{n+1} = y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

Where the slopes are:
*   $k_1 = h f(x_n, y_n)$
*   $k_2 = h f(x_n + \frac{h}{2}, y_n + \frac{k_1}{2})$
*   $k_3 = h f(x_n + \frac{h}{2}, y_n + \frac{k_2}{2})$
*   $k_4 = h f(x_n + h, y_n + k_3)$

## Algorithm
1.  **Define function**: $f(x, y)$.
2.  **Input**: Initial $x_0, y_0$, target $x_n$, and step size $h$.
3.  **Iterate**: Until current $x$ reaches target $x_n$:
    *   Calculate $k_1, k_2, k_3, k_4$.
    *   Calculate weighted average $k = (k_1 + 2k_2 + 2k_3 + k_4) / 6$.
    *   Update $y = y + k$.
    *   Update $x = x + h$.
4.  **Output**: Final value of $y$.

## Pseudocode
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
