#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

#define INFINITY_ROOTS -1
#define EPS 1e-7

int eq_solver(double a, double b, double c, double* x1, double* x2);
bool is_equal(double a, double b);
void test_eq_solver();

int main()
{
    double a, b, c, x1 = 0, x2 = 0;

    cout << "Quadratic equation solver" << endl << endl;
    cout << "Input coefficients of equation a*x^2 + b*x + c = 0" << endl;
    cout << endl << "a = ";
    cin >> a;
    cout << endl << "b = ";
    cin >> b;
    cout << endl << "c = ";
    cin >> c;
    cout << endl;

    int num_roots = eq_solver(a, b, c, &x1, &x2);

    if(num_roots == INFINITY_ROOTS)
        cout << "Infinity roots" << endl;
    else if(num_roots == 0)
        cout << "No roots" << endl;
    else if(num_roots == 1)
        cout << "One root: x = " << x1 << endl;
    else
        cout << "Two roots: x1 = " << x1 << ", x2 = " << x2 << endl;
    return 0;
}


/*!
    Solves the quadratic equation ax^2 + bx + c = 0
    /param a  [in]  a-coefficient
    /param b  [in]  b-coefficient
    /param c  [in]  c-coefficient
    /param x1 [out] Pointer to the first root
    /param x2 [out] Pointer to the second root

    /return Number of roots

    /note If number of roots is infinite, returns INFINITY_ROOTS
*/

int eq_solver(double a, double b, double c, double* x1, double* x2)
{
    assert(isfinite(a));
    assert(isfinite(b));
    assert(isfinite(c));
    assert(x1 != NULL);
    assert(x2 != NULL);
    assert(x1 != x2);

    if(is_equal(a, 0) && is_equal(b, 0)){
        if(is_equal(c, 0))
            return INFINITY_ROOTS;
        else
            return 0;
    }
    if(is_equal(a, 0)){
        *x1 = -(c/b);
        return 1;
    }


    double D = b*b - 4*a*c;
    if(D < 0){
        return 0;
    }
    else if(is_equal(D, 0)){
        *x1 = -(b / (2*a));
        return 1;
    }
    else{
        *x1 = (-b - sqrt(D)) / (2*a);
        *x2 = (-b + sqrt(D)) / (2*a);
        return 2;
    }
}


/*!
    Checks two float or double numbers on equality

    /param [in] a First number
    /param [in] b Second number

    /return True if a is equal to b or false otherwise

    /note Comparison is done with precision EPS
*/

bool is_equal(const double a, const double b)
{
    assert(isfinite(a));
    assert(isfinite(b));

    if(a > b - EPS && a < b + EPS)
        return true;
    else
        return false;
}


void test_eq_solver()
{
    double x1 = 0, x2 = 0;

    assert(eq_solver(0, 0, 0, &x1, &x2) == INFINITY_ROOTS);
    assert(eq_solver(0, 0, 0.000001, &x1, &x2) == 0);
    assert(eq_solver(0, 0, 0.00000001, &x1, &x2) == INFINITY_ROOTS);
    assert(eq_solver(0, -3, 1, &x1, &x2) == 1 && is_equal(x1, (double)1/3));
    assert(eq_solver(1, 1, 1, &x1, &x2) == 0);
    assert(eq_solver(2, 4, 2, &x1, &x2) == 1 && is_equal(x1, -1));
    assert(eq_solver(1, 5, 4, &x1, &x2) == 2 && is_equal(x1, -4) && is_equal(x2, -1));
}
