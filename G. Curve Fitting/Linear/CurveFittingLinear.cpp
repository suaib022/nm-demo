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
