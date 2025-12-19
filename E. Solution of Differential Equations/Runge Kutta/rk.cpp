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
