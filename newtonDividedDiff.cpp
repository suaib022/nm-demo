#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

double calculateProterm(int i, double value, vector<double>& x) {
    double pro = 1;
    for (int j = 0; j < i; j++) {
        pro = pro * (value - x[j]);
    }
    return pro;
}

int main()
{
    int n;
    cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n);
    vector<vector<double>> y(n, vector<double>(n));

    cout << "Enter data points (x and y):" << endl;
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i][0];
    }

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++) {
            y[j][i] = (y[j + 1][i - 1] - y[j][i - 1]) / (x[i + j] - x[j]);
        }
    }

    cout << endl << "Divided Difference Table:" << endl;
    for (int i = 0; i < n; i++) {
        cout << setw(10) << fixed << setprecision(2) << x[i];
        for (int j = 0; j < n - i; j++) {
            cout << setw(10) << fixed << setprecision(4) << y[i][j];
        }
        cout << endl;
    }

    double value;
    cout << endl << "Enter the value to interpolate: ";
    cin >> value;

    double sum = y[0][0];
    for (int i = 1; i < n; i++) {
        sum = sum + (calculateProterm(i, value, x) * y[0][i]);
    }

    cout << endl << "Value at " << value << " is " << fixed << setprecision(4) << sum << endl;
    double initial_ans = sum;

    char choice;
    cout << endl << "Use additional data point? (y/n): ";
    cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        n++;
        x.resize(n);
        y.resize(n, vector<double>(n));

        cout << "Enter the new data point (x and y): ";
        cin >> x[n - 1] >> y[n - 1][0];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < n - i; j++) {
                y[j][i] = (y[j + 1][i - 1] - y[j][i - 1]) / (x[i + j] - x[j]);
            }
        }

        cout << endl << "Updated Divided Difference Table:" << endl;
        for (int i = 0; i < n; i++) {
            cout << setw(10) << fixed << setprecision(2) << x[i];
            for (int j = 0; j < n - i; j++) {
                cout << setw(10) << fixed << setprecision(4) << y[i][j];
            }
            cout << endl;
        }

        sum = y[0][0];
        for (int i = 1; i < n; i++) {
            sum = sum + (calculateProterm(i, value, x) * y[0][i]);
        }

        cout << endl << "New Value at " << value << " is " << fixed << setprecision(4) << sum << endl;
        cout << "Error: " << abs(sum - initial_ans) << endl;
    }

    return 0;
}

/*
Sample Input/Output:

Enter number of data points: 4
Enter data points (x and y):
5 12
6 13
9 14
11 16

Divided Difference Table:
      5.00   12.0000    1.0000   -0.1667    0.0500
      6.00   13.0000    0.3333    0.1333
      9.00   14.0000    1.0000
     11.00   16.0000

Enter the value to interpolate: 10

Value at 10.0000 is 14.6667

Use additional data point? (y/n): y
Enter the new data point (x and y): 12 18

Updated Divided Difference Table:
      5.00   12.0000    1.0000   -0.1667    0.0500   -0.0024
      6.00   13.0000    0.3333    0.1333    0.0333
      9.00   14.0000    1.0000    0.3333
     11.00   16.0000    2.0000
     12.00   18.0000

New Value at 10.0000 is 14.7143
Error: 0.0476
*/
