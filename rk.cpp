#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Differential equation: dy/dx = x + y
double f(double x, double y) {
    return x + y;
}

int main() {
    double x0, y0, xn, h;

    cout << "Enter initial values (x0, y0): ";
    cin >> x0 >> y0;
    
    cout << "Enter target x (xn): ";
    cin >> xn;
    
    cout << "Enter step size (h): ";
    cin >> h;

    cout << endl << fixed << setprecision(4);
    
    int n = (int)((xn - x0) / h);
    double y = y0;
    double x = x0;

    for (int i = 0; i <= n; i++) {
        cout << "x" << i << " = " << x << ", y" << i << " = " << y << endl;

        if (i == n) break; 

        double k1 = h * f(x, y);
        double k2 = h * f(x + h / 2.0, y + k1 / 2.0);
        double k3 = h * f(x + h / 2.0, y + k2 / 2.0);
        double k4 = h * f(x + h, y + k3);

        double k = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

        y = y + k;
        x = x + h;
    }
    
    cout << endl << "Value of y at x = " << xn << " is " << y << endl;

    return 0;
}

/*
Sample Input/Output:

Enter initial values (x0, y0): 0 1
Enter target x (xn): 0.2
Enter step size (h): 0.1

x0 = 0.0000, y0 = 1.0000
x1 = 0.1000, y1 = 1.1103
x2 = 0.2000, y2 = 1.2428

Value of y at x = 0.2000 is 1.2428
*/