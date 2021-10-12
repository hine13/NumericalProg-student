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

## Integration

### rectangularRect()

Integral using rectangular Method.

**Parameters**

* **x\[\]**: input data vector **x** in 1D-array .
* **y\[\]**: input data vector **y** in 1D-array.
* **m**: The number associated with the number of for statements.

**Example code**

```c
double IntegrateRect(double x[], double y[], int m) {
	int N = m - 1;
	double I = 0;
	for (int i = 0; i < N; i++)
		I += y[i] * (x[i + 1] - x[i]);

	return I;
}
```



### integral\(\)

Integral using trapezoidal Method.

**Parameters**

* **x\[\]**: input data vector **x** in 1D-array .
* **y\[\]**: input data vector **y** in 1D-array.
* **m**: The number associated with the number of for statements.

**Example code**

```c
double trapz(double x[], double y[], int m) {
	int N = m - 1;
	double I = 0;
	for (int i = 0; i < N; i++) {
		I += (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
	}
	return I * 0.5;
}
```



Integral using Simpson 1/3 Method.

```text
double integral(double func(const double _x), double a, double b, int n);
```

**Parameters**

* **func(const double x)** : A small function that goes into a function.
* **x\[\]**: input data vector **x** in 1D-array .
* **y\[\]**: input data vector **y** in 1D-array.
* **dydx\[\]**: output vector **dydx** in 1D-array.
* **m**:  length **x** and **y**.

**code**

```c
double integral(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double I = func(a) + 4 * func(b - h) + func(b);
	for (int i = 1; i < n - 2; i += 2)
	{
		double xi = a + i * h;
		I += 4 * func(xi) + 2 * func(xi + h);
	}
	return I * h / 3;
}
```



Integral using Simpson 3/8 Method.

**Parameters**

* **func**: Function **func** is defined.
* **a** is starting point of x.
* **b** is ending point of x.
* **n** is the length between **a** and **b**

**Example code**

```c
double integral38(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double I = func(a) + 3 * func(b - h) + 3 * func(b - 2 * h) + func(b);
	for (int i = 1; i < n - 4; i += 3)
	{
		double xi = a + i * h;
		I += 3 * func(xi) + 3 * func(xi + h) + 2 * func(xi + 2 * h);
	}
	return I * h * 3 / 8;
}
```



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

**Example code**

```c
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	double* t;
	t = (double*)malloc(sizeof(double));
	y[0] = 0;
	t[0] = 0;
	int N = (tf - t0) / h;
	printf("odeEU\n");
	for (int i = 0; i < N; i++) {
		double slope = myfunc(t[i], y[i]);
		t[i + 1] = t[i] + h;
		y[i + 1] = y[i] + slope * h;
		printf("%f\n", y[i]);
	}

	free(t);
}
```

## odeEM\(\)

Solve the 1st-order ODE using Euler's modified Method.

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: Solution of ODE in structure 1D-array form.
* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.

**Example code**

```c
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h) {
	double* t;
	t = (double*)malloc(sizeof(double));
	y[0] = 0;
	t[0] = 0;
	int N = (tf - t0) / h;
	printf("odeEM\n");
	for (int i = 0; i < N; i++) {
		double slope1 = myfunc(t[i], y[i]);
		t[i + 1] = t[i] + h;
		y[i + 1] = y[i] + slope1 * h;
		double slope2 = myfunc(t[i + 1], y[i + 1]);
		y[i + 1] = y[i] + (myfunc(t[i], y[i]) + myfunc(t[i + 1], y[i + 1])) * h / 2.0;
		printf("%f\n", y[i]);
	}

	free(t);

}
```

#  odeRK2\(\)

Solve the 1st-order ODE using Runge-Kutta second order Method.

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: Solution of ODE in structure 1D-array form.
* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.

**Example code**

```c
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	double* t;
	double c1 = 0.5;
	double c2 = 0.5;
	double a2 = 1.0;
	double b21 = 1.0;
	t = (double*)malloc(sizeof(double));
	y[0] = 0;
	t[0] = 0;
	int N = (tf - t0) / h;
	printf("odeRK2\n");
	for (int i = 0; i < N; i++) {
		double K1 = myfunc(t[i], y[i]);
		t[i + 1] = t[i] + h;
		double K2 = myfunc(t[i] + a2 * h, y[i] + (K1 * h) * b21);
		y[i + 1] = y[i] + ((c1 * K1 + c2 * K2) * h);
		printf("%f\n", y[i]);
	}



	free(t);

}
```

# odeRK3

Solve the 1st-order ODE using Runge-Kutta third order Method.

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: Solution of ODE in structure 1D-array form.
* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.

**Example code**

```c
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	double* t;
	double c1 = 1.0 / 6.0;
	double c2 = 4.0 / 6.0;
	double c3 = 1.0 / 6.0;
	double a2 = 0.5;
	double a3 = 1.0;
	double b21 = 0.5;
	double b31 = -1.0;
	double b32 = 2.0;
	t = (double*)malloc(sizeof(double));
	y[0] = 0;
	t[0] = 0;
	int N = (tf - t0) / h;
	printf("odeRK3\n");
	for (int i = 0; i < N; i++) {
		double K1 = myfunc(t[i], y[i]);
		t[i + 1] = t[i] + h;
		double K2 = myfunc(t[i] + h * a2, y[i] + (K1 * h) * b21);
		double K3 = myfunc(t[i] + a3 * h, y[i] + b31 * K1 * h + b32 * K2 * h);
		y[i + 1] = y[i] + ((K1 + 4 * K2 + K3) * h) / 6.0;
		printf("%f\n", y[i]);
	}



	free(t);

}
```

 # ode 

Also, create a function that calls different ODE method

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: Solution of ODE in structure 1D-array form.
* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **method**: choose method

**Example code**

```c
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method) {
	switch (method)
	{
	case 0:
		odeEM(myfunc, y, t0, tf, h);
		break;
	case 1:
		odeEU(myfunc, y, t0, tf, h);
		break;
	case 2:
		odeRK2(myfunc, y, t0, tf, h);
		break;
	case 3:
		odeRK3(myfunc, y, t0, tf, h);
		break;
	}

}
```

 # ode part2

Also, create a function that calls different ODE method

**Parameters**

* **myfunc1,2**: Function **func** is defined.
* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y1_init**: y[0]
* **y2_init**: z[0]

**Example code**

```c
void sys2RK2(double myfunc1(const double t, const double y, double z), double myfunc2(const double t, const double y, double z), double t0, double tf, double h, double y1_init, double y2_init) {
	double N = (tf - t0) / h;
	double* t;
	t = (double*)malloc(sizeof(double)*N);
	double* y;
	y = (double*)malloc(sizeof(double) * N);
	double* z;
	z = (double*)malloc(sizeof(double) * N);
	y[0] = y1_init;
	z[0] = y2_init;
	t[0] = t0;
	
	printf("t[0]  %f  %f  %f\n", t[0], y[0], z[0]);

	for (int i = 0; i < N+1; i++) {
		t[i + 1] = t[i] + h;

		double ky1 = myfunc1(t[i], y[i], z[i]);
		double kz1 = myfunc2(t[i], y[i], z[i]);
		double ky2 = myfunc1(t[i] + h, y[i] + ky1 * h, z[i] + kz1 * h);
		double kz2 = myfunc2(t[i] + h, y[i] + ky1 * h, z[i] + kz1 * h);

		y[i + 1] = y[i] + h * (ky1 + ky2) / 2.0;
		z[i + 1] = z[i] + h * (kz1 + kz2) / 2.0;
		
		printf("t[%d]  %f  %f  %f\n", i, t[i], y[i], z[i]);

	}
	free(t);
	free(y);
	free(z);
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
