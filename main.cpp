#include <iostream>
#include <cmath>
#include <stdio.h>
#define RADIUS 5   
#define WIN_W 1000
#define WIN_H 1000

const int dim = 3;
const double RAD_TO_ARCSECONDS = 206264.8;
const double GMBH = 6.743015e-11 * 1.988416e30 * 4.15e6;	// G * black hole mass
const double RBH = 8.53e3 * 206265 * 1.495978707e11;
const double C = 299792458.0;
const double Y = 31558149.0;			// sideric period
const double T = Y * 20;		
const double step = Y / 1000.0;
double const min_h = Y / 1000000.0;
const int N = 1;


void f_GR(double* X, double* m, double* dX)
{
	for (int i = 0; i < dim; i++) {
		dX[i] = X[i + dim * N];
	}
	double R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	double Vsqrd = X[3] * X[3] + X[4] * X[4] + X[5] * X[5];
	double VR = X[0] * X[3] + X[1] * X[4] + X[2] * X[5];
	for (int i = 0; i < dim; i++) {
		dX[i + 3] = -GMBH / (C * C * R * R * R) * (X[i] * (C * C - 4 * GMBH / R + Vsqrd) - 4 * X[i + 3] * VR);
	}
}

void Ralston3(double* X, double* Xbuf, double* k1, double* k2, double* k3, double* m, double h)
{
	f_GR(X, m, k1);

	for (int i = 0; i < dim * 2 * N; i++) {
		Xbuf[i] = X[i] + h * 0.5 * k1[i];
	}
	f_GR(Xbuf, m, k2);

	for (int i = 0; i < dim * 2 * N; i++) {
		Xbuf[i] = X[i] + h * 0.75 * k2[i];
	}
	f_GR(Xbuf, m, k3);

	for (int i = 0; i < dim * 2 * N; i++) {
		X[i] += h * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]);
	}
}

double distance(double* v1, double* v2) {
	return std::sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) + pow(v1[2] - v2[2], 2));
}

double scalar(double* a, double* b) {
	double res = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return res;
}

void calculateEquatorialCoordinates(double x, double y, double z, double* result) {
	// ��������� ������ ����������� (alpha) � ��������
	double alpha = atan2(y, x); // atan2 ��������� ����� x � y

	// ��������� ��������� (delta) � ��������
	double r = sqrt(x * x + y * y + z * z); // ���������� �� ������ ���������
	double delta = asin(z / r); // ��������� = arcsin(z / r)

	// ����������� ���������� �� ������ � ������� ����
	result[0] = alpha * RAD_TO_ARCSECONDS; // ������ ����������� � �������� ����
	result[1] = delta * RAD_TO_ARCSECONDS; // ��������� � �������� ����
}

void calculateCartesianCoordinates(double alpha_arcseconds, double delta_arcseconds, double r, double* result) {
	// ��������� ������ ����������� � ��������� �� ���������� � �������
	double alpha = alpha_arcseconds / RAD_TO_ARCSECONDS;
	double delta = delta_arcseconds / RAD_TO_ARCSECONDS;

	// ��������� ��������� ����������
	result[0] = r * cos(delta) * cos(alpha); // x
	result[1] = r * cos(delta) * sin(alpha); // y
	result[2] = r * sin(delta);              // z
}


int main()
{
	bool exactFlag = false;
	double h = min_h;
	double t = 0;		
	double m[N];
	double X[dim * N * 2];
	double Xbuf[dim * N * 2];
	double k1[dim * N * 2];
	double k2[dim * N * 2];
	double k3[dim * N * 2];
	double k4[dim * N * 2];
	double EqCoords[2] = {0, 0};
	double EqCoordsBH[2] = {959100.7014, -104377.4682};		// R.A., Decl in arcsec
	double BH[3] = {0, 0, 0};

	calculateCartesianCoordinates(EqCoordsBH[0], EqCoordsBH[1], RBH, BH);


	X[0] = -13946410030006.707;		// S2
	X[1] = 2625500133928.8228;		// T_P = 2002.32
	X[2] = 12102295968349.498;

	X[3] = 2689.819254116237e3;
	X[4] = 6762.698965156799e3;
	X[5] = 1632.5708144507482e3;


	FILE* fp_S2 = fopen("S2.txt", "w");
	FILE* fp_S2_Eq = fopen("S2_Equatorial.txt", "w");
	t = 0.32 * Y;
	while (t < T)
	{
		for (int i = 0; i < step / h; i++, t += h) {
			Ralston3(X, Xbuf, k1, k2, k3, m, h);
		}
		calculateEquatorialCoordinates(X[0] + BH[0], X[1] + BH[1], X[2] + BH[2], EqCoords);

		fprintf(fp_S2_Eq, "%.3f %.4f %.4f\n", t / Y + 2002, EqCoordsBH[0] - EqCoords[0] - 360 * 3600, EqCoordsBH[1] - EqCoords[1]);
		fprintf(fp_S2, "%.3f %.16le %.16le %.16le\n", t / Y + 2002, X[0], X[1], X[2]);
	}
	fclose(fp_S2);
	fclose(fp_S2_Eq);

	X[0] = 3310508584282.186;		// S38
	X[1] = 30953372574500.473;		// T_P = 2003.30
	X[2] = 2347171645064.925;

	X[3] = 5374.5976777876895e3;
	X[4] = -670.1138144615753e3;
	X[5] = 1256.6745270510417e3;

	FILE* fp_S38 = fopen("S38.txt", "w");
	FILE* fp_S38_Eq = fopen("S38_Equatorial.txt", "w");
	t = 1.30 * Y;
	while (t < T)
	{
		for (int i = 0; i < step / h; i++, t += h) {
			Ralston3(X, Xbuf, k1, k2, k3, m, h);
		}
		calculateEquatorialCoordinates(X[0] + BH[0], X[1] + BH[1], X[2] + BH[2], EqCoords);

		fprintf(fp_S38_Eq, "%.3f %.4f %.4f\n", t / Y + 2002, EqCoordsBH[0] - EqCoords[0] - 360 * 3600, EqCoordsBH[1] - EqCoords[1]);
		fprintf(fp_S38, "%.3f %.16le %.16le %.16le\n", t / Y + 2002, X[0], X[1], X[2]);
	}
	fclose(fp_S38);
	fclose(fp_S38_Eq);

	X[0] = 30495208238659.08;		// S55
	X[1] = -5657388258855.3955;		// T_P = 2009.31
	X[2] = 15610192040343.342;

	X[3] = 262.34894244566254e3;
	X[4] = -4657.414935485226e3;
	X[5] = -2200.4335446310647e3;

	FILE* fp_S55 = fopen("S55.txt", "w");
	FILE* fp_S55_Eq = fopen("S55_Equatorial.txt", "w");
	t = 7.31 * Y;
	while (t < T)
	{
		for (int i = 0; i < step / h; i++, t += h) {
			Ralston3(X, Xbuf, k1, k2, k3, m, h);
		}
		calculateEquatorialCoordinates(X[0] + BH[0], X[1] + BH[1], X[2] + BH[2], EqCoords);

		fprintf(fp_S55_Eq, "%.3f %.4f %.4f\n", t / Y + 2002, EqCoordsBH[0] - EqCoords[0] - 360 * 3600, EqCoordsBH[1] - EqCoords[1]);
		fprintf(fp_S55, "%.3f %.16le %.16le %.16le\n", t / Y + 2002, X[0], X[1], X[2]);
	}
	fclose(fp_S55);
	fclose(fp_S55_Eq);
}