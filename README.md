# Spline-Interpolation
This python scripts interpolate the Runge-Function and a curve function with natural kubic Splines.
Please check the results at the bottom.


We want to interpolate a function with a given dataset (xi,fi), i = 1,...,n

Firstly, the python script builds a polynom of the form si(x) = c1,i+c2,i(x-xi)+c3,i(x-xi)^2+c4,i(x-xi)^3 for each intervall [xi,xi+1], i = 1,...,n-1
Secondly, the spline function is built, which contains si of every intervall.


The spline function accepts 4 arguments.
The dataset xi and fi, a precomputed S = (S1, ..., Sn), the derivatives at the points xi, which is relevant to calculate the coefficients and of course the x value for which we want to evaluate the y value.

For the kubic splines, we only need to determine the second derivative at x1 and xn to conclude all Si by solving a linear equation system.
Since we deal with natural splines, we use the fact that s1''(x1) = 0 and s2''(xn) = 0 to be able to solve the system.
We precompute S with the funtion ```solve_S(xi,fi,0,0)```.
If different conditions are necessary for s1''(x1) and s2''(xn), simply set them as second and third argument here.

After precomputing 

```S = solve_S(xi,fi,0,0)```,
we can now evaluate the spline function simply with

```spline(xi,fi,S,x)```,
where x ∈ [x1,xn] is the value we want to evaluate.


![Spline auf Runge mit n=7](https://github.com/nokitoino/Spline-Interpolation/blob/main/myplot.png "n=7 datasets")
  n=7
![Spline auf Kurve mit n=15](https://github.com/nokitoino/Spline-Interpolation/blob/main/myplot1.png "n=15 datasets")
  n=15
![Spline auf Kurve mit n=9](https://github.com/nokitoino/Spline-Interpolation/blob/main/myplot2.png "n=9 datasets")
  n=9
