# Newton-Raphson method
## used to get the root approximation of the
## real-valued function (wikipedia.org)

import math

def f(x):
    """
    define the y as a function of x
    let say f(x) = e^x - 5x^2
    """

    return math.exp(x) - (5 * (x ** 2))

def diff_f(x):
    """
    define the differential or derivative of y
    as a function of x
    in case of y = f(x) = e^x - 5x^2, dy/dx = e^x - 10x
    """

    return math.exp(x) - (10 * x)

def newton_raphson(x, fx, dfx):
    """
    Newton-Raphson method
    will be get the new x according to the slope
    of the differencial function
    """
    
    return x - (fx / dfx)

iter = 0
epsilon = 1E-15
x0 = 1

while (True):
    dfx = (f(x0 + 1E-10) - f(x0)) / 1E-10
    x1 = newton_raphson(x0, f(x0), dfx)

    if (abs(x1 - x0) < epsilon):
        break
    
    iter += 1
    x0 = x1

print("x =", round(x0, 6), "with", iter, "iterations")

def bisection(a, b, eps):
    """
    bisection method
    
    input:
        a: int | float
        b: int | float
        eps: float

    output:
        c: int: float
    """
    a, b = a, b

    while (abs(a - b) >= eps):
        c = (a + b) / 2
        fa, fc = f(a), f(c)

        if (fa * fc < 0):
            b = c

        else:
            a = c

    return c