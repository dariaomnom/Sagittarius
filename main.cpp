#include <iostream>
#include <cmath>
#include <stdio.h>
#define RADIUS 5   
#define WIN_W 1000
#define WIN_H 1000

const int dim = 3;
const double PI = 3.14159265358979323846;
const double RAD_TO_ARCSECONDS = 206264.8;
const double GMBH = 6.743015e-11 * 1.988416e30 * 4.15e6;	// G * black hole mass
const double RBH = 8.53e3 * 206265 * 1.495978707e11;
const double C = 299792458.0;
const double Y = 31558149.0;			// sideric period
const double T = Y * 2024;
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

void multiplyMatrixVector(const double matrix[9], double vector[3]) {
	double result[3];
	for (int i = 0; i < 3; ++i) {
		result[i] = 0;
		for (int j = 0; j < 3; ++j) {
			result[i] += matrix[i * 3 + j] * vector[j];
		}
	}
	vector[0] = result[0];
	vector[1] = result[1];
	vector[2] = result[2];
}

void multiplyMatrixMatrix(double matrix1[9], double matrix2[9]) {
	double result[9];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			result[i * 3 + j] = 0;
			for (int k = 0; k < 3; ++k) {
				result[i * 3 + j] += matrix1[i * 3 + k] * matrix2[k * 3 + j];
			}
		}
	}
	for (int i = 0; i < 9; ++i) {
		matrix2[i] = result[i];
	}
}

void transpose(double matrix[9]) {
	double buffer;
	buffer = matrix[1];
	matrix[1] = matrix[3];
	matrix[3] = buffer;

	buffer = matrix[2];
	matrix[2] = matrix[6];
	matrix[6] = buffer;

	buffer = matrix[5];
	matrix[5] = matrix[7];
	matrix[7] = buffer;
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
	double EqCoords[2] = { 0, 0 };				// ����� ��� ��������� ������
	double EqCoordsBH[2] = { 959100.7014, -104377.4682 };	// R.A., Decl in arcsec
	double BH[3] = { 0, 0, 0 }; 				// ��������� ���������� ��
	//double RotateMatrix1[9] = {					
	//	0, 0, 1,
	//	0, 1, 0,
	//	-1, 0, 0
	//};
	//double RotateMatrix2[9] = {
	//	1, 0, 0,
	//	0, std::cos(-959100.7014 / RAD_TO_ARCSECONDS), -std::sin(-959100.7014 / RAD_TO_ARCSECONDS),
	//	0, std::sin(-959100.7014 / RAD_TO_ARCSECONDS), std::cos(-959100.7014 / RAD_TO_ARCSECONDS)
	//};
	//double RotateMatrix3[9] = {
	//	std::cos(-104377.4682 / RAD_TO_ARCSECONDS), 0, std::sin(-104377.4682 / RAD_TO_ARCSECONDS),
	//	0, 1, 0,
	//	-std::sin(-104377.4682 / RAD_TO_ARCSECONDS), 0, std::cos(-104377.4682 / RAD_TO_ARCSECONDS)
	//};

	double RotateMatrix0[9] = {
	0, 1, 0,
	1, 0, 0,
	0, 0, 1
	};
	double RotateMatrix1[9] = {					
		std::cos(EqCoordsBH[0] / RAD_TO_ARCSECONDS), -std::sin(EqCoordsBH[0] / RAD_TO_ARCSECONDS), 0,
		std::sin(EqCoordsBH[0] / RAD_TO_ARCSECONDS), std::cos(EqCoordsBH[0] / RAD_TO_ARCSECONDS), 0,
		0, 0, 1
	};
	double RotateMatrix2[9] = {
		std::cos(PI / 2 - EqCoordsBH[1] / RAD_TO_ARCSECONDS), 0, std::sin(PI / 2 - EqCoordsBH[1] / RAD_TO_ARCSECONDS),
		0, 1, 0,
		-std::sin(PI / 2 - EqCoordsBH[1] / RAD_TO_ARCSECONDS), 0, std::cos(PI / 2 - EqCoordsBH[1] / RAD_TO_ARCSECONDS)
	};
	double RotateMatrix3[9] = {
	0, 1, 0,
	-1, 0, 0,
	0, 0, 1
	};

	multiplyMatrixMatrix(RotateMatrix1, RotateMatrix0);
	multiplyMatrixMatrix(RotateMatrix2, RotateMatrix0);
	multiplyMatrixMatrix(RotateMatrix3, RotateMatrix0);
	//transpose(RotateMatrix1);

	calculateCartesianCoordinates(EqCoordsBH[0], EqCoordsBH[1], RBH, BH);   // ������ ���������� ��������� �� ��������������, ������������ � ������ BH


	X[0] = -13946410030007.033;		// ���������� S2
	X[1] = 2625500133928.027;		// ������ ����������� ���������� T_P = 2002.32
	X[2] = 12102295968349.293;
	multiplyMatrixVector(RotateMatrix0, X);

	X[3] = 2689819.2541161072;		// �������� S2
	X[4] = 6762698.965156819;
	X[5] = 1632570.8144508821;
	multiplyMatrixVector(RotateMatrix0, X + 3);


	FILE* fp_S2 = fopen("S2.txt", "w");
	FILE* fp_S2_Eq = fopen("S2_Equatorial.txt", "w");
	t = 2002.32 * Y;			// ������ ������� ������ �������
	while (t < T)
	{
		for (int i = 0; i < step / h; i++, t += h) {
			Ralston3(X, Xbuf, k1, k2, k3, m, h);
		}
		calculateEquatorialCoordinates(X[0] + BH[0], X[1] + BH[1], X[2] + BH[2], EqCoords); 				// ���������� �������������� ���������. X - ���������� ������ � ������� ������������ ��

		// if ((2010.0 - t/Y) < 0.001){
			fprintf(fp_S2_Eq, "%.3f %.4f %.4f\n", t / Y, EqCoordsBH[0] - EqCoords[0] - 360 * 3600, EqCoordsBH[1] - EqCoords[1]);		// � ����: ������� ����� ���. ������������ �� � ������
			fprintf(fp_S2, "%.3f %.16le %.16le %.16le\n", t / Y, X[0], X[1], X[2]);						// � ����: ����, ��������� ���������� ������������ ��
		// }
	}
	fclose(fp_S2);
	fclose(fp_S2_Eq);

	X[0] = 3447941175082.65;		// S38
	X[1] = 32238372168557.91;		// T_P = 2003.30
	X[2] = 2444612226178.7485;
	multiplyMatrixVector(RotateMatrix0, X);

	X[3] = 5364230.710015342;
	X[4] = -668821.2436060756;
	X[5] = 1254250.5494617736;
	multiplyMatrixVector(RotateMatrix0, X + 3);

	FILE* fp_S38 = fopen("S38.txt", "w");
	FILE* fp_S38_Eq = fopen("S38_Equatorial.txt", "w");
	t = 2003.30 * Y;
	while (t < T)
	{
		for (int i = 0; i < step / h; i++, t += h) {
			Ralston3(X, Xbuf, k1, k2, k3, m, h);
		}
		calculateEquatorialCoordinates(X[0] + BH[0], X[1] + BH[1], X[2] + BH[2], EqCoords);

		// if ((2010.0 - t/Y) < 0.001){
			fprintf(fp_S38_Eq, "%.3f %.4f %.4f\n", t / Y, EqCoordsBH[0] - EqCoords[0] - 360 * 3600, EqCoordsBH[1] - EqCoords[1]);
			fprintf(fp_S38, "%.3f %.16le %.16le %.16le\n", t / Y, X[0], X[1], X[2]);
		// }
	}
	fclose(fp_S38);
	fclose(fp_S38_Eq);

	X[0] = 31761187579458.13;		// S55
	X[1] = -5892249309894.6;		// T_P = 2009.31
	X[2] = 16258234200748.242;
	multiplyMatrixVector(RotateMatrix0, X);

	X[3] = 261842.9021437616;
	X[4] = -4648431.32900279;
	X[5] = -2196189.166724118;
	multiplyMatrixVector(RotateMatrix0, X + 3);

	FILE* fp_S55 = fopen("S55.txt", "w");
	FILE* fp_S55_Eq = fopen("S55_Equatorial.txt", "w");
	t = 2009.31 * Y;
	while (t < T)
	{
		for (int i = 0; i < step / h; i++, t += h) {
			Ralston3(X, Xbuf, k1, k2, k3, m, h);
		}
		calculateEquatorialCoordinates(X[0] + BH[0], X[1] + BH[1], X[2] + BH[2], EqCoords);

		// if ((2010.0 - t/Y) < 0.0001) {
			fprintf(fp_S55_Eq, "%.3f %.4f %.4f\n", t / Y, EqCoordsBH[0] - EqCoords[0] - 360 * 3600, EqCoordsBH[1] - EqCoords[1]);
			fprintf(fp_S55, "%.3f %.16le %.16le %.16le\n", t / Y, X[0], X[1], X[2]);
		// }
	}
	fclose(fp_S55);
	fclose(fp_S55_Eq);
}
