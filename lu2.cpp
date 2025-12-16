#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
    char choice;
    do {
        int n;
        cout << "Enter number of equations: ";
        cin >> n;

        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);

        cout << "Enter the augmented matrix:" << endl;
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) cin >> a[i][j];
            cin >> b[i];
        }

        vector<vector<double>> l(n, vector<double>(n, 0));
        vector<vector<double>> u(n, vector<double>(n, 0));

        for(int i=0; i<n; i++) {
            for(int k=i; k<n; k++) {
                double sum = 0;
                for(int j=0; j<i; j++) sum += l[i][j] * u[j][k];
                u[i][k] = a[i][k] - sum;
            }

            for(int k=i; k<n; k++) {
                if(i == k) l[i][i] = 1;
                else {
                    double sum = 0;
                    for(int j=0; j<i; j++) sum += l[k][j] * u[j][i];
                    
                    if(abs(u[i][i]) > 1e-9)
                        l[k][i] = (a[k][i] - sum) / u[i][i];
                }
            }
        }

        cout << "\nThe lower and upper triangular matrix" << endl;
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++)
                cout << fixed << setprecision(3) << l[i][j] << "  ";
            cout << endl;
        }
        cout << endl;

        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++)
                cout << fixed << setprecision(3) << u[i][j] << "  ";
            cout << endl;
        }

        bool unique = true;
        for(int i=0; i<n; i++) {
            if(abs(u[i][i]) < 1e-9) {
                unique = false;
                break;
            }
        }

        vector<double> y(n);
        for(int i=0; i<n; i++) {
            double sum = 0;
            for(int j=0; j<i; j++) sum += l[i][j] * y[j];
            y[i] = b[i] - sum;
        }

        if(unique) {
            vector<double> x(n);
            for(int i=n-1; i>=0; i--) {
                double sum = 0;
                for(int j=i+1; j<n; j++) sum += u[i][j] * x[j];
                x[i] = (y[i] - sum) / u[i][i];
            }

            cout << "\nThe solution for each variable rounded to 3 decimal places" << endl;
            for(int i=0; i<n; i++)
                cout << "X" << i+1 << " = " << fixed << setprecision(3) << x[i] << endl;
            
            cout << "\nThe system has unique solution" << endl;
        } else {
            bool no_solution = false;
            for(int i=n-1; i>=0; i--) {
                if(abs(u[i][i]) < 1e-9 && abs(y[i]) > 1e-9) {
                    no_solution = true;
                    break;
                }
            }
            if(no_solution) cout << "\nThe system has no solution" << endl;
            else cout << "\nThe system has infinite solution" << endl;
        }

        cout << "Solve another system? (y/n): ";
        cin >> choice;
    } while(choice == 'y' || choice == 'Y');

    return 0;
}

/*
Sample Input/Output:

Enter number of equations: 5
Enter the augmented matrix:
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15

The lower and upper triangular matrix
1.000  0.000  0.000  0.000  0.000
0.500  1.000  0.000  0.000  0.000
1.500  0.200  1.000  0.000  0.000
1.000  0.000  0.800  1.000  0.000
0.500  -0.600  0.800  1.714  1.000

2.000  1.000  -1.000  3.000  2.000
0.000  2.500  2.500  -2.500  0.000
0.000  0.000  5.000  -3.000  -5.000
0.000  0.000  0.000  1.400  3.000
0.000  0.000  0.000  0.000  1.857

The solution for each variable rounded to 3 decimal places
X1 = 5.154
X2 = -1.000
X3 = 2.262
X4 = -0.138
X5 = 1.185

The system has unique solution
Solve another system? (y/n): n
*/