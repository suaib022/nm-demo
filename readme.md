# Numerical Methods

## Summary
*   [LU Decomposition Method](#lu-decomposition-method)
*   [Newton Forward Interpolation](#newton-forward-interpolation)
*   [Newton Backward Interpolation](#newton-backward-interpolation)
*   [Newton Divided Difference Interpolation](#newton-divided-difference-interpolation)
*   [Runge-Kutta (RK) Method](#runge-kutta-rk-method)

---

<a id="lu-decomposition-method"></a>
# LU Decomposition Method

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

## Pseudocode
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

---

<a id="newton-divided-difference-interpolation"></a>
# Newton Divided Difference Interpolation

## Theory
Newton's Divided Difference Interpolation is a method for constructing an interpolating polynomial for a given set of data points where the interval between data points is **not necessarily equal**.

The formula is given by:
$$f(x) = f(x_0) + (x-x_0)f[x_0, x_1] + (x-x_0)(x-x_1)f[x_0, x_1, x_2] + \dots + (x-x_0)\dots(x-x_{n-1})f[x_0, \dots, x_n]$$

Where divided differences are defined recursively:
$$f[x_i, x_j] = \frac{f(x_j) - f(x_i)}{x_j - x_i}$$

## Algorithm
1.  **Input**: Read ($n$) data points ($x, y$) and the value to interpolate ($value$).
2.  **Table**: Construct the divided difference table using the recursive formula.
3.  **Compute**: Initialize $sum = y[0][0]$.
    *   Iterate $i$ from 1 to $n-1$.
    *   Accumulate product term: $\prod_{j=0}^{i-1} (value - x_j)$
    *   Add ($product \times y[0][i]$) to sum.
4.  **Output**: Display the interpolated value.

## Pseudocode
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
