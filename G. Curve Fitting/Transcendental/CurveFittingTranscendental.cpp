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
