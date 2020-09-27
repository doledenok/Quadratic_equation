#include <iostream>
#include <cmath>
#include <assert.h>

/**
    Used for description of infinity number of roots in equation
*/
#define INFINITY_ROOTS -1

/**
    Used for comparison of two float numbers
*/
#define EPS 1e-7


/**
    Solves the quadratic equation ax^2 + bx + c = 0

    \param a  [in]  a-coefficient
    \param b  [in]  b-coefficient
    \param c  [in]  c-coefficient
    \param x1 [out] Pointer to the first root
    \param x2 [out] Pointer to the second root

    \return Number of roots

    \note If number of roots is infinite, returns INFINITY_ROOTS
*/
int eq_solver(double a, double b, double c, double* x1, double* x2);


/**
    Checks two float or double numbers on equality

    \param a [in] First number
    \param b [in] Second number

    \return True if a is equal to b or false otherwise

    \note Comparison is done with precision EPS
*/
bool is_equal(double a, double b);


/**
    Checks work of function "eq_solver" with asserts
*/
void test_eq_solver();

int main()
{
    double a, b, c, x1 = 0, x2 = 0;

    std::cout << "Quadratic equation solver\n\n";
    std::cout << "Input coefficients of equation a*x^2 + b*x + c = 0\n";
    std::cout << "\na = ";
    std::cin >> a;
    std::cout << "\nb = ";
    std::cin >> b;
    std::cout << "\nc = ";
    std::cin >> c;
    std::cout << std::endl;

    int num_roots = eq_solver(a, b, c, &x1, &x2);

    if(num_roots == INFINITY_ROOTS)
        std::cout << "Infinity roots\n";
    else if(num_roots == 0)
        std::cout << "No roots\n";
    else if(num_roots == 1)
        std::cout << "One root: x = " << x1 << std::endl;
    else
        std::cout << "Two roots: x1 = " << x1 << ", x2 = " << x2 << std::endl;
    return 0;
}


int eq_solver(double a, double b, double c, double* x1, double* x2)
{
    assert(std::isfinite(a));
    assert(std::isfinite(b));
    assert(std::isfinite(c));
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


bool is_equal(const double a, const double b)
{
    assert(std::isfinite(a));
    assert(std::isfinite(b));

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
