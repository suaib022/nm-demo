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
