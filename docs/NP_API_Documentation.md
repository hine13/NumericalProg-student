---
description: Numerical Method example
---


# Reference of Numerical Programming API

`#include "myNM.h"`



# SinTaylor

Create a simple function that returns the output of sine x.

```c
double sinTaylor(double _x)
```

**parameters**

* **N_max**, **epsilon** : related with stop condition
* **S_N_prev** : previous S_N
* **N** : The number associated with the number of for statements.

**Example code**

```c
double sinTaylor(double _x)
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;

	do {
		N++;
		S_N_prev = S_N;
		S_N = 0;
		for (int k = 0; k < N; k++)
			S_N += pow(-1, k) * pow(_x, 2 * k + 1) / factorial(2 * k + 1);

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	return S_N;
}
```

```c
double sinTaylor2(double _x)
```

**parameters**

* **N_max**, **epsilon** : related with stop condition
* **S_N_prev** : previous S_N
* **N** : The number associated with the number of for statements.

**Example code**

```c
double sinTaylor2(double _x)
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;

	do {
		N++;
		S_N_prev = S_N;
		S_N = 0;
		for (int k = 0; k < N; k++)
		{
			//  (-1)^n			
			int sign_part = 1;
			for (int i = 1; i <= k; i++)   //sign_part *= -1;
				sign_part *= -1;

			//  (x)^n			
			double pow_part = 1;
			for (int i = 1; i <= 2 * k + 1; i++)  //pow_part *= _x * _x;
				pow_part *= _x;

			// Factorial
			double fac_part = 1;
			for (int i = 1; i <= 2 * k + 1; i++)
				fac_part *= i;

			S_N += sign_part * pow_part / fac_part;
		}

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	return S_N;
}
```



## Non-Linear Solver

### newtonRaphson\(\)

Solves the non-linear problem using Newton-Raphson method

```text
double newtonRaphson(double x0, double tol);
```

**Parameters**

* **x0:**  initial value.
* **tol**:  tolerance error

**Example code**

```text
double tol = 0.00001;
double x0 = 3;
double NR_result;

NR_result = newtonRaphson(x0, tol);
```

## Linear Solver

### gaussElim\(\)

solves for vector **x** from Ax=b, a linear system problem Using Gauss Elimination

```cpp
void gaussElim(Matrix _A, Matrix _B, Matrix* _U, Matrix* _B_out);
```

**Parameters**

* **A**: Matrix **A** in structure Matrix form. Should be \(nxn\) square.
* **B**: vector **b** in structure Matrix form. Should be \(nx1\)
* **U**: Matrix **U** in structure Matrix form. Should be \(nxn\) square.
* **B\_out**: vector **B\_out** in structure Matrix form. Should be \(nx1\)

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix vecd = zeros(vecb.rows, vecb.cols);

gaussElim(matA, vecb, matU, vecd);
```

### inv\(\)

Find the inverse Matrix.

```text
void inv(Matrix _A, Matrix _Ainv);
```

**Parameters**

* **A**: Matrix **A** in structure Matrix form. Should be \(nxn\) square.
* **Ainv**: Matrix **Ainv** in structure Matrix form. Should be \(nxn\) square.

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix matAinv = zeros(matA.rows, matA.cols);

inv(matA, matAinv);
```

### 

## Numerical Differentiation

### gradient1D\(\)

Solve for numerical gradient \(dy/dt\) from a 1D-array form.

```text
void gradient1D(double x[], double y[], double dydx[], int m);
```

**Parameters**

* **x\[\]**: input data vector **x** in 1D-array .
* **y\[\]**: input data vector **y** in 1D-array.
* **dydx\[\]**: output vector **dydx** in 1D-array.
* **m**:  length **x** and **y**.

**code**

```cpp
void gradient1D(double x[], double y[], double dydx[], int m) {
	double h = x[1] - x[0];

	if (sizeof(x) != sizeof(y)) {
		printf("ERROR: length of x and y must be equal\n");
		return;
	}

	// three-Point FWD  O(h)
	// Modify to Three-Point FWD  O(h^2)
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);

	// Two-Point Central  O(h^2)
	for (int i = 1; i < m - 1; i++) {
		dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);
	}

	// Two-Point BWD  O(h)
	// Modify to Three-Point BWD  O(h^2)
	dydx[m - 1] = (3 * y[m - 1] - 4 * y[m - 2] + y[m - 3]) / (2 * h);
```

### gradientfunc\(\)

How can you use ‘myFunc()’ in the function of ‘gradientFunc()’?

**Parameters**

* **func(const double x)** : A small function that goes into a function.
* **x\[\]**: input data vector **x** in 1D-array .
* **y\[\]**: input data vector **y** in 1D-array.
* **dydx\[\]**: output vector **dydx** in 1D-array.
* **m**:  length **x** and **y**.

**code**

```c
void gradientFunc(double func(const double x), double x[], double dydx[], int m) {
	double* y;

	y = (double*)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++) {
		y[i] = func(x[i]);
	}

	//printVec(y, m);
	gradient1D(x, y, dydx, m);

	free(y);
}
```



Solve for numerical g

## Integration

### integral\(\)

Integral using Simpson 1/3 Method.

```text
double integral(double func(const double _x), double a, double b, int n);
```

**Parameters**

* **func**: Function **func** is defined.
* **a** is starting point of x.
* **b** is ending point of x.
* **n** is the length between **a** and **b**

**Example code**

```text
double I_simpson13 = integral(myFunc, -1, 1, 12);

double myFunc(const double _x) {
	return sqrt(1 - (_x * _x));
}
```

## ODE-IVP

### odeEU\(\)

Solve the 1st-order ODE using Euler's Explicit Method.

```text
void odeEU(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: Solution of ODE in structure 1D-array form.
* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y0** is initial value of **y\[\]**.

**Example code**

```text
double a = 0;
double b = 0.1;
double h = 0.001;
double y_EU[200] = { 0 };
double v0 = 0;

odeEU(odeFunc_rc, y_EU, a, b, h, v0);

double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double omega = 2 * PI * f;
	return  -T * v + T * Vm * cos(omega * t);
}
```

-------------------------------------------------------------------------------------------------------

## Class or Header name

### Function Name

```text

```

**Parameters**

* p1
* p2

**Example code**

```text

```
