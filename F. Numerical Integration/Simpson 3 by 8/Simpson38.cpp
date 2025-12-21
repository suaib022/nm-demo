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
