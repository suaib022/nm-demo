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

vector<double> rangess(double step)
{
    double dmax = 0.0;
    for (int j = 0; j <= degree; j++)
        dmax = max(dmax, fabs(coef[j] / coef[0]));

    double low = -(1 + dmax);
    double high = (1 + dmax);
    vector<double> res;

    while (low < high)
    {
        if (f(low) * f(low + step) < 0)
        {
            if (fabs(f(low)) < fabs(f(low + step)))
                res.push_back(low);
            else
                res.push_back(low + step);
        }
        low += step;
    }
    return res;
}

double newton_raphson(double initial_guess, double tol)
{
    double x0 = initial_guess;
    double fx = f(x0);
    if (fabs(fx) < 1e-12)
        return x0;

    double dfx = df(x0);
    if (fabs(dfx) <= 1e-12)
        return NAN;

    double x1 = x0 - fx / dfx;

    while (fabs(x1 - x0) >= tol)
    {
        x0 = x1;
        fx = f(x0);
        if (fabs(fx) < 1e-12)
            return x0;
        dfx = df(x0);
        if (fabs(dfx) <= 1e-12)
            return NAN;
        x1 = x0 - fx / dfx;
    }
    return x1;
}

int main()
{
    ifstream fin("inputnewtraph.txt");
    ofstream fout("outputnewtraph.txt");

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

        vector<double> initial_guesses = rangess(0.45);
        double tol = 1e-6;
        vector<double> roots;

        for (double guess : initial_guesses)
        {
            double root = newton_raphson(guess, tol);
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
