/**********************************************************

Code written in C to calculate the stress fields of an
elliptical inclusion. We assume plane strain conditions.
The eigen strain is generic: non-dilatational with shear
components. This implementation is based on

Two-dimensional elastic inclusion problems,
M A Jaswon and R D Bhargava,
Mathematical Proceedings of the Cambridge Philosophical Society,
Vol. 53, Issue 7, July 1961, pp. 669-680.

Wherever equation numbers are used in the comments, they refer
to this paper cited above.

Author: M P Gururajan
Date: 22, April 2013
***********************************************************/
/***********************************************************
Updated by Nitin Davessar
Made changes to the equations where ever needed.
Added comments and references for all the equations used.
Date: 26 April, 2013
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

int main(void)
{

	double kappa;
	double nu;
	double c1, c2, c3, c4;
	double mu, mu1, lambda, lambda1, e, e1;
	double A[3][3], B[3];
	double Omega11, Omega22, Omega12;
	double eps[4];
	double a, b, c, d;
	double complex x, xbar;
	double complex f;
	double psidp;
	double complex denom;
	double tmp;
	double epsilon11, epsilon22, epsilon12;
	int i, j;
	double X, Y;
	double sigxx, sigyy, sigxy;
	double sigxxs, sigyys, sigxys;
	double chi, eta;
	double aplusb;
	double C1;
	double complex C2;
	double T1, T2, T3a, T3b;
	double T4, T5, T6, T7;
	double complex T8, T9;
	double complex T10, T11;
	double theta;
	double complex K1, K2;
	double sig11, sig22, sig12;
	double sig11s, sig22s, sig12s;
	double complex eq21a, eq21b;
	FILE *fp1;
	FILE *input;

	// Input parameters from input file
	input = fopen("inputs", "r");
	// Read the input parameters
	char buf[100];
	fscanf(input, "%s ", buf);
	nu = atof(buf); /* Poisson's ratio */
	fscanf(input, "%s ", buf);
	a = atof(buf); /* Major axis */
	fscanf(input, "%s ", buf);
	b = atof(buf); /* Minor axis */
	fscanf(input, "%s ", buf);
	e = atof(buf); /* Young's Modulus of Matrix */
	fscanf(input, "%s ", buf);
	e1 = atof(buf); /* Young's Modulus of Inclusion */
	fscanf(input, "%s ", buf);
	Omega11 = atof(buf); /* Omega11 */
	fscanf(input, "%s ", buf);
	Omega22 = atof(buf); /* Omega22 */
	fscanf(input, "%s ", buf);
	Omega12 = atof(buf); /* Omega12 */

	mu = e / (2 * (1 + nu));					   /*  The shear modulus of the matrix phase */
	lambda = e * nu / ((1 + nu) * (1 - 2 * nu));   /* The other Lame's constant for the matrix phase */
	mu1 = e1 / (2 * (1 + nu));					   /* The shear modulus of the precipitate phase */
	lambda1 = e1 * nu / ((1 + nu) * (1 - 2 * nu)); /*The other Lame's constant for the precipitate phase */
	kappa = 3.0 - 4.0 * nu;						   /* Plane strain condition */

	/* In this part of the code, we calculate the equivalent eigenstrain */
	/* This part is probably taken from Mura, Micromechanics of defects in solids */

	c1 = b / ((1.0 + kappa) * (a + b)) * (2.0 + kappa + ((a - b) / (a + b)));
	c2 = b / ((1.0 + kappa) * (a + b)) * (2.0 - kappa - ((a - b) / (a + b)));
	c3 = a / ((1.0 + kappa) * (a + b)) * (2.0 - kappa + ((a - b) / (a + b)));
	c4 = a / ((1.0 + kappa) * (a + b)) * (2.0 + kappa - ((a - b) / (a + b)));

	A[1][1] = 2.0 * mu * (c1 - 1.0) + lambda * (c1 + c3 - 1.0) - 2.0 * mu1 * c1 - lambda1 * (c1 + c3);
	A[1][2] = 2.0 * mu * c2 + lambda * (c2 + c4 - 1.0) - 2.0 * mu1 * c2 - lambda1 * (c2 + c4);
	A[2][1] = 2.0 * mu * c3 + lambda * (c1 + c3 - 1.0) - 2.0 * mu1 * c3 - lambda1 * (c1 + c3);
	A[2][2] = 2.0 * mu * (c4 - 1.0) + lambda * (c2 + c4 - 1.0) - 2.0 * mu1 * c4 - lambda1 * (c2 + c4);

	B[1] = -2.0 * mu1 * Omega11 - lambda1 * (Omega11 + Omega22);
	B[2] = -2.0 * mu1 * Omega22 - lambda1 * (Omega11 + Omega22);

	eps[1] = (A[2][2] * B[1] - A[1][2] * B[2]) / (A[1][1] * A[2][2] - A[1][2] * A[2][1]);
	eps[2] = (A[2][1] * B[1] - A[1][1] * B[2]) / (A[1][2] * A[2][1] - A[1][1] * A[2][2]);

	eps[3] = mu1 * Omega12 /
			 (mu - (mu - mu1) * (1.0 - ((4.0 * a * b) / ((1.0 + kappa) * (a + b) * (a + b)))));

	epsilon11 = eps[1];
	epsilon22 = eps[2];
	epsilon12 = eps[3];

	printf("Equivalent eigen strain: %le %le %le\n", eps[1], eps[2], eps[3]);

	/* Calculation of stress fields */
	/* Note that in Jaswon and Bhargava, the stress fields are given */
	/* for principal eigenstrains and pure shear eigenstrains */
	/* Hence, in what follows below, we also do these two calculations */
	/* separately, and add the same for the most generic inclusion case */

	fp1 = fopen("analytical.dat", "w");
	c = sqrt(a * a - b * b);
	if (a == b)
		printf("c is zero. You might get some NaNs\n");
	aplusb = a + b;

	for (i = 0; i < 750; ++i)
	{
		chi = (double)i * 0.01;
		for (j = 0; j < 22; ++j)
		{
			eta = j * 0.314;
			X = c * cosh(chi) * cos(eta);
			Y = c * sinh(chi) * sin(eta);
			K1 = exp(chi) * (cos(eta) + _Complex_I * sin(eta)) - exp(-chi) * (cos(eta) - _Complex_I * sin(eta));
			K2 = exp(chi) * (cos(eta) - _Complex_I * sin(eta)) - exp(-chi) * (cos(eta) + _Complex_I * sin(eta));
			theta = 0.5 * acos(creal(K1 / K2));
			if ((X * X / (a * a) + Y * Y / (b * b)) <= 1.0)
			{
				/* Principal eigenstrain part: Equation 14  */
				/* Note: sig11, sig22 and sig12 are, respectively, */
				/* P11, P22 and P12 of Equation 14 */
				sig12 = 0.0;
				sig11 = (b + 2.0 * a) * epsilon11 / aplusb + b * epsilon22 / aplusb;
				sig11 = -4.0 * mu * a * sig11 / ((1.0 + kappa) * aplusb);
				sig22 = a * epsilon11 / aplusb + (a + 2.0 * b) * epsilon22 / aplusb;
				sig22 = -4.0 * mu * b * sig22 / ((1.0 + kappa) * aplusb);
				/* Shear eigenstrain part: Equation 22 */
				/* Here sig11s, sig22s and sig12s are P11, P22 and P12 respectivley*/
				sig22s = 0.0;
				sig11s = 0.0;
				sig12s = -8 * mu * epsilon12 * a * b / ((kappa + 1.0) * (aplusb) * (aplusb));
				/* Total stress = principal + shear*/
				sig11 = sig11 + sig11s;
				sig22 = sig22 + sig22s;
				sig12 = sig12 + sig12s;
			}
			else
			{
				/* Principal eigenstrain part: Equation 13 */
				/* Note: sigxx, sigyy and sigxy are, respectively, */
				/* p(chi)(chi), p(eta)(eta) and p(chi)(eta) */
				C1 = 1.0 - (sinh(2.0 * chi) / (cosh(2.0 * chi) - cos(2.0 * eta)));
				C1 = 8.0 * mu * a * b * (epsilon11 - epsilon22) * C1 / ((1.0 + kappa) * c * c);
				T1 = cosh(2.0 * chi) - ((a * a + b * b) / (c * c));
				T2 = sinh(2.0 * chi);
				T3a = 2.0 * a * a / (c * c);
				T3b = 2.0 * b * b / (c * c);
				T4 = exp(chi);
				T5 = cos(eta);
				T6 = sin(eta);
				T7 = exp(-chi);
				T8 = T4 * (T5 + _Complex_I * T6) + T7 * (T5 - _Complex_I * T6);
				T8 = T8 / (T4 * (T5 + _Complex_I * T6) - T7 * (T5 - _Complex_I * T6));
				T9 = 1.0 - exp(-2.0 * chi) * cos(2.0 * eta) +
					 _Complex_I * exp(-2.0 * chi) * sin(2.0 * eta);
				T10 = T1 * T8 - T2 + T3a * T9;
				T11 = -T1 * T8 + T2 - T3b * T9;
				C2 = 8.0 * mu * a * b / ((1.0 + kappa) * c * c * (cosh(2.0 * chi) - cos(2.0 * eta)));
				C2 = C2 * (T10 * epsilon11 + T11 * epsilon22);
				sigxy = 0.5 * cimag(C2);
				sigxx = 0.5 * (C1 - creal(C2));
				sigyy = 0.5 * (C1 + creal(C2));
				/* We transform confocal co-ordinates back to Cartestian - Equation 6 */
				/* Note: sig11, sig22, sig12 are p11, p22 and p12 respectively*/
				sig22 = 0.5 * (sigxx + sigyy + (sigyy - sigxx) * cos(2.0 * theta) - 2.0 * sigxy * sin(2.0 * theta));
				sig11 = 0.5 * (sigxx + sigyy - (sigyy - sigxx) * cos(2.0 * theta) + 2.0 * sigxy * sin(2.0 * theta));
				sig12 = 0.5 * (2.0 * sigxy * cos(2.0 * theta) - (sigyy - sigxx) * sin(2.0 * theta));
				/* Shear eigenstrain part: Equation 21 */
				/* Here, sigxxs, sigyys and sigxys are */
				/* p(chi)(chi), p(eta)(eta) and p(chi)(eta) respectively*/
				C1 = sin(2.0 * eta) / (cosh(2.0 * chi) - cos(2.0 * eta));
				C2 = -8. * mu * epsilon12 * 2. * a * b / ((1. + kappa) * c * c);
				T1 = cosh(2. * chi) - ((a * a + b * b) / (c * c));
				T8 = T1 / tanh(chi + _Complex_I * eta);
				T9 = -sinh(2. * chi) + ((a * a + b * b) / (c * c)) * (1. - exp(-2. * (chi + i * eta)));
				T10 = T8 + T9;
				T11 = T10 / (cosh(2.0 * chi) - cos(2.0 * eta));
				T11 = _Complex_I * -C2 * T11;
				eq21a = C2 * C1;
				eq21b = __real__(T11);
				sigyys = 0.5 * (eq21a + eq21b);
				sigxxs = 0.5 * (eq21a - eq21b);
				sigxys = 0.5 * __imag__(T11);
				/* We transform confocal co-ordinates back to Cartestian - Equation 6 */
				/* Note: sig11s, sig22s, sig12s are p11, p22 and p12 respectively*/
				sig22s = 0.5 * (sigxxs + sigyys + (sigyys - sigxxs) * cos(2.0 * theta) - 2.0 * sigxys * sin(2.0 * theta));
				sig11s = 0.5 * (sigxxs + sigyys - (sigyys - sigxxs) * cos(2.0 * theta) + 2.0 * sigxys * sin(2.0 * theta));
				sig12s = 0.5 * (2.0 * sigxys * cos(2.0 * theta) - (sigyys - sigxxs) * sin(2.0 * theta));
				/* Total stress = principal + shear */
				sig11 = sig11 + sig11s;
				sig22 = sig22 + sig22s;
				sig12 = sig12 + sig12s;
			}
			if (j == 0)
				fprintf(fp1, "%le %le %le %le %le\n", X, Y, sig11, sig22, sig12);
		}
	}
	fclose(fp1);
	return 0;
}
