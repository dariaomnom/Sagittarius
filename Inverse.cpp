#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <unistd.h>
// #include <Windows.h>
// #include <unistd>

#define ITER
#define DEBUG
const int stars_num = 3; // S2, S38, S55
const int states_num = 6 * 3; // (3 coords + 3 velocities) * 3 stars
const int params_num = states_num + 1; // states + mass BH
const int observation_num = 33 + 29 + 26; // alpha, delta for each observation of three stars
const double RAD_TO_ARCSECONDS = 206264.8;
const double C = 299792458.0;
const double start_time = 1995.0;
const double end_time = 2016.0;
const double Y = 31558149.0;
const double STEP = Y / 10000.0;
const std::vector<double> BH_cart = { -1.4388065248297531e+19, -2.2977127354968464e+20, -1.2758089026580264e+20 };
const std::vector<double> BH_eq = { 959100.7014, -104377.4682 };


// FILE* fp_S2 = fopen("S2_debug.txt", "w");
// 				FILE* fp_S38 = fopen("S38_debug.txt", "w");
// 				FILE* fp_S55 = fopen("S55_debug.txt", "w");

std::vector<double> F(std::vector<double>& condition) {
	std::vector<double> result(params_num + states_num * params_num, 0.0),
		R(stars_num),
		Vsqrd(stars_num),
		VxR(stars_num),
		dFdP(states_num / 2, 0.0),
		dFdx(states_num * states_num, 0.0);

	double GM = condition[params_num - 1];

	R[0] = sqrt(condition[0] * condition[0] + condition[1] * condition[1] + condition[2] * condition[2]);
	R[1] = sqrt(condition[3] * condition[3] + condition[4] * condition[4] + condition[5] * condition[5]);
	R[2] = sqrt(condition[6] * condition[6] + condition[7] * condition[7] + condition[8] * condition[8]);

	Vsqrd[0] = condition[9] * condition[9] + condition[10] * condition[10] + condition[11] * condition[11];
	Vsqrd[1] = condition[12] * condition[12] + condition[13] * condition[13] + condition[14] * condition[14];
	Vsqrd[2] = condition[15] * condition[15] + condition[16] * condition[16] + condition[17] * condition[17];

	VxR[0] = condition[0] * condition[9] + condition[1] * condition[10] + condition[2] * condition[11];
	VxR[1] = condition[3] * condition[12] + condition[4] * condition[13] + condition[5] * condition[14];
	VxR[2] = condition[6] * condition[15] + condition[7] * condition[16] + condition[8] * condition[17];

	// Derivatives of coordinates
	for (int i = 0; i < states_num / 2; i++) {
		result[i] = condition[i + states_num / 2];
	}

	// Derivatives of velocities
	for (int i = 0; i < states_num / 2; i++) {
		result[i + states_num / 2] = -GM / (C * C * R[i / 3] * R[i / 3] * R[i / 3]) \
			* (condition[i] * (C * C - 4 * GM / R[i / 3] + Vsqrd[i / 3]) - 4 * condition[i + states_num / 2] * VxR[i / 3]);
	}

	// Unit matrix 9x9 in right upper block of dF/dx
	for (int i = 0; i < states_num / 2; i++) {
		dFdx[states_num * i + states_num / 2 + i] = 1;
	}

	// left lower block of dF/dx
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				int ik = 9 + i + 3 * k;
				int jk = j + 3 * k;
				if (ik - 9 == jk)
					dFdx[states_num * ik + jk] = -GM * (pow(R[k], 2) - 3 * pow(condition[jk], 2)) / pow(R[k], 5);
				else
					dFdx[states_num * ik + jk] = GM * 3 * condition[ik - 9] * condition[jk] / pow(R[k], 5);
			}
		}
	}

	// result product of martices
	for (int i = 0; i < states_num; i++) { // (dF/dP) + (dF/dx)(dx/dP)
		for (int j = 0; j < params_num; j++) {

			if (j == params_num - 1 && i >= states_num / 2)
				result[i * params_num + j + params_num] = -condition[i - states_num / 2] / pow(R[(i - states_num / 2) / 3], 3.);

			for (int k = 0; k < states_num; k++) {
				result[i * params_num + j + params_num] += dFdx[i * states_num + k] * condition[k * params_num + j + params_num];
			}
		}
	}
	return result;
}


void RK4_step(std::vector<double>& condition, std::vector<std::vector<double> >& bufers, double h)
{
	std::vector<double> conditionBuf = bufers[0], k1 = bufers[1], k2 = bufers[2], k3 = bufers[3];
	k1 = F(condition);
	for (int i = 0; i < params_num + states_num * params_num; i++) {
		conditionBuf[i] = condition[i] + 0.5 * h * k1[i];
	}

	k2 = F(conditionBuf);
	for (int i = 0; i < params_num + states_num * params_num; i++) {
		conditionBuf[i] = condition[i] + 0.75 * h * k2[i];
	}

	k3 = F(conditionBuf);
	for (int i = 0; i < params_num + states_num * params_num; i++) {
		condition[i] += h * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]);
	}
}


double g(std::vector<double>& condition, int observ_counter, int k, int func_num) { // testing

	double x = condition[3 * k] + BH_cart[0];
	double y = condition[3 * k + 1] + BH_cart[1];
	double z = condition[3 * k + 2] + BH_cart[2];
	double r = sqrt(x * x + y * y + z * z);

	if (func_num == 0) {
		double alpha = atan2(y, x);
		return alpha * RAD_TO_ARCSECONDS - BH_eq[0];
	}
	if (func_num == 1) {
		double delta = asin(z / r);
		return delta * RAD_TO_ARCSECONDS - BH_eq[1];
	}
	return 0;
}


std::vector<double> dgdx(std::vector<double>& condition, int observ_counter, int k, int func_num) { // testing
	std::vector<double> result(states_num, 0.0);

	double x = condition[3 * k] + BH_cart[0];
	double y = condition[3 * k + 1] + BH_cart[1];
	double z = condition[3 * k + 2] + BH_cart[2];
	double r = sqrt(x * x + y * y + z * z);

	if (func_num == 0) {
		result[k * 3] = -y / (x * x + y * y) * RAD_TO_ARCSECONDS;
		result[k * 3 + 1] = x / (x * x + y * y) * RAD_TO_ARCSECONDS;
	}
	if (func_num == 1) {
		result[k * 3] = -z * x / (pow(r, 3.0) * pow(1 - pow(z / r, 2.0), 1.0 / 2)) * RAD_TO_ARCSECONDS;
		result[k * 3 + 1] = -z * y / (pow(r, 3.0) * pow(1 - pow(z / r, 2.0), 1.0 / 2)) * RAD_TO_ARCSECONDS;
		result[k * 3 + 2] = (x * x + y * y) / (pow(r, 3.0) * pow(1 - pow(z / r, 2.0), 1.0 / 2)) * RAD_TO_ARCSECONDS;
	}

	return result;
}


// Функция для выполнения LUP-разложения
void LUPDecompose(std::vector<double>& A, std::vector<int>& P, int n) {
	for (int i = 0; i < n; ++i) P[i] = i;

	for (int k = 0; k < n; ++k) {
		// Поиск максимального элемента в текущем столбце
		int maxRow = k;
		for (int i = k + 1; i < n; ++i) {
			if (abs(A[i * n + k]) > abs(A[maxRow * n + k])) {
				maxRow = i;
			}
		}

		// Меняем строки местами в перестановочном векторе
		std::swap(P[k], P[maxRow]);

		// Меняем строки матрицы
		for (int j = 0; j < n; ++j) {
			std::swap(A[k * n + j], A[maxRow * n + j]);
		}

		// Выполняем преобразования для разложения
		for (int i = k + 1; i < n; ++i) {
			A[i * n + k] /= A[k * n + k];
			for (int j = k + 1; j < n; ++j) {
				A[i * n + j] -= A[i * n + k] * A[k * n + j];
			}
		}
	}
}

// Прямая подстановка для решения Ly = Pb
void forwardSubstitution(const std::vector<double>& A, std::vector<double>& y, const std::vector<double>& Pb, int n) {
	for (int i = 0; i < n; ++i) {
		y[i] = Pb[i];
		for (int j = 0; j < i; ++j) {
			y[i] -= A[i * n + j] * y[j];
		}
	}
}

// Обратная подстановка для решения Ux = y
void backSubstitution(const std::vector<double>& A, std::vector<double>& x, const std::vector<double>& y, int n) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = y[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= A[i * n + j] * x[j];
		}
		x[i] /= A[i * n + i];
	}
}

// Основная функция для решения СЛАУ методом LUP-разложения
std::vector<double> solveLUP(const std::vector<double>& A, const std::vector<double>& B, int n) {
	std::vector<double> LU = A;       // Рабочая копия матрицы A
	std::vector<int> P(n);            // Перестановочный вектор
	std::vector<double> Pb(n), y(n), x(n);

	// Разложение LUP
	LUPDecompose(LU, P, n);

	// Преобразуем вектор B с учетом перестановок
	for (int i = 0; i < n; ++i) {
		Pb[i] = B[P[i]];
	}

	// Решаем Ly = Pb
	forwardSubstitution(LU, y, Pb, n);

	// Решаем Ux = y
	backSubstitution(LU, x, y, n);

	return x;
}



std::vector<double> Gauss_Newton(std::vector<double>& params, std::vector<double>& observations, std::vector<std::vector<double> >& bufers) {
	std::vector<double> condition(params_num + states_num * params_num, 0.0),
		b(params_num, 0.0),                                 // AT * Wsqrd * r(beta)
		x(params_num, 0.0),                                 // AT * Wsqrd * r(beta)
		B(params_num * params_num, 0.0),                    // (AT * W * A)
		dg(states_num, 0.0),
		Ak(params_num, 0.0);

				FILE* fp_S2 = fopen("S2_debug.txt", "w");
				FILE* fp_S38 = fopen("S38_debug.txt", "w");
				FILE* fp_S55 = fopen("S55_debug.txt", "w");

	for (int i = 0; i < params_num; i++) condition[i] = params[i];
	for (int i = 0; i < states_num; i++) {
		condition[i * params_num + i + params_num] = 1.0;
	}
	double t = start_time * Y;
	double step = STEP;
	int observ_counter = 0;
	while (t < end_time * Y and observ_counter < observation_num) {

		if (observations[observ_counter * 6] * Y - t < step) {
			step = observations[observ_counter * 6] * Y - t;
			// step = 0.0;
			RK4_step(condition, bufers, step);

				// fprintf(fp_S2, "%.3f %.16le %.16le %.16le\n", t / Y, condition[0], condition[1], condition[2]);
				// fprintf(fp_S38, "%.3f %.16le %.16le %.16le\n", t / Y, condition[3], condition[4], condition[5]);
				// fprintf(fp_S55, "%.3f %.16le %.16le %.16le\n", t / Y, condition[6], condition[7], condition[8]);

			t += step;
#ifdef ITER 
			std::cout << "\n";
			std::cout << "Observation number = " << observ_counter << std::endl;
			std::cout << "Time = " << observations[observ_counter * 6] << std::endl;
			std::cout << "step = " << step << std::endl;
#endif

#ifdef DEBUG 
			std::cout << "\n  CONDITION  \n" << std::endl;
			for (int i = 0; i < condition.size(); i++) {
				std::cout << std::setprecision(1) << std::setw(7) << condition[i] << " ";
				if ((i + 1) % 19 == 0) printf("\n");
			}
			printf("\n");
			// sleep(1);
#endif


			// RA
			dg = dgdx(condition, observ_counter, observations[observ_counter * 6 + 5], 0);

#ifdef DEBUG 
			std::cout << "\n  D G  \n" << std::endl;
			for (int i = 0; i < dg.size(); i++) {
				std::cout << std::setprecision(1) << std::setw(7) << dg[i] << " ";
				if ((i + 1) % 19 == 0) printf("\n");
			}
			printf("\n");
			// sleep(1);
#endif

			for (int i = 0; i < params_num; i++) {	// line of A matrix -(dg/dx * dx/dP)_k
				Ak[i] = 0.0;
				for (int j = 0; j < states_num; j++) {
					Ak[i] -= dg[j] * condition[j * params_num + i + params_num];
				}
			}

#ifdef DEBUG 
			std::cout << "\n  Ak \n" << std::endl;
			for (int i = 0; i < Ak.size(); i++) {
				std::cout << std::setprecision(0) << std::setw(7) << Ak[i] << " ";
				if ((i + 1) % 19 == 0) printf("\n");
			}
			printf("\n");
			// sleep(1);
#endif

			for (int i = 0; i < params_num; i++) {	// A * Wsqrd * r
				b[i] += Ak[i] / observations[observ_counter * 6 + 3] * (observations[observ_counter * 6 + 1] - g(condition, observ_counter, observations[observ_counter * 6 + 5], 0));
			}
			// b = At * sqrt(w) * r(beta)
			// At - один столбец k
			// r - одна невязка r[i]

			for (int i = 0; i < params_num; i++) {	// AT * W * A
				for (int j = 0; j < params_num; j++) {
					B[i * params_num + j] += Ak[i] * Ak[j] / pow(observations[observ_counter * 6 + 3], 2.0);
					if (Ak[i] == 0 || Ak[j] == 0)  B[i * params_num + j] = 0.0;
#ifdef EBUG 
					std::cout << "\nB[i * params_num + j] = " << B[i * params_num + j] << "= " << Ak[i] << " * " << Ak[j] << " / " << pow(observations[observ_counter * 6 + 3], 2.0) << std::endl;
					// sleep(1);
#endif
				}
			}

			// Decl
			dg = dgdx(condition, observ_counter, observations[observ_counter * 6 + 5], 1);
			for (int i = 0; i < params_num; i++) {	// line of A matrix -(dg/dx * dx/dP)_k
				Ak[i] = 0.0;
				for (int j = 0; j < states_num; j++) {
					Ak[i] -= dg[j] * condition[j * params_num + i + params_num];
				}
			}

			for (int i = 0; i < params_num; i++) {	// A * Wsqrd * r
				b[i] += Ak[i] / observations[observ_counter * 6 + 4] * (observations[observ_counter * 6 + 2] - g(condition, observ_counter, observations[observ_counter * 6 + 5], 1));
			}

			for (int i = 0; i < params_num; i++) {	// AT * W * A
				for (int j = 0; j < params_num; j++) {
					B[i * params_num + j] += Ak[i] * Ak[j] / pow(observations[observ_counter * 6 + 4], 2.0);
				}
			}

			observ_counter++;
			step = STEP;
			/*#ifdef DEBUG
				for (int i = 0; i < params_num; i++) {
					for (int j = 0; j < params_num; j++) {
						printf("%.3f ", B[i * params_num + j]);
					}
					printf("\n");
				}
			}*/
			/*#ifdef DEBUG
				for (int i = 0; i < params_num; i++) {
					printf("%.3f ", b[i]);
				}
				printf("\n");
			}*/
		}
		else {
			RK4_step(condition, bufers, step);
				fprintf(fp_S2, "%.3f %.16le %.16le %.16le\n", t / Y, condition[0], condition[1], condition[2]);
				fprintf(fp_S38, "%.3f %.16le %.16le %.16le\n", t / Y, condition[3], condition[4], condition[5]);
				fprintf(fp_S55, "%.3f %.16le %.16le %.16le\n", t / Y, condition[6], condition[7], condition[8]);
			t += step;
		}

	}
#ifdef DEBUG  // матрица B
	// printf("\n   B MATRIX\n");
	// for (int i = 0; i < 19 * 19; i++) {
	// 	printf("%.0f ", B[i]);
	// 	if ((i + 1) % 19 == 0) printf("\n");
	// }
	// printf("\n");

			std::cout << "\n  B MATRIX  \n" << std::endl;
			for (int i = 0; i < B.size(); i++) {
				std::cout << std::setprecision(1) << std::setw(7) << B[i] << " ";
				if ((i + 1) % 19 == 0) printf("\n");
			}
			printf("\n");
			// sleep(1);
#endif
#ifdef DEBUG 
	// printf("\n   b MATRIX\n");
	// for (int i = 0; i < params_num; i++) {
	// 	printf("%.3f ", b[i]);
	// }
	// printf("\n\n");
			std::cout << "\n  b MATRIX  \n" << std::endl;
			for (int i = 0; i < b.size(); i++) {
				std::cout << std::setprecision(1) << std::setw(7) << b[i] << " ";
				if ((i + 1) % 19 == 0) printf("\n");
			}
			printf("\n");

#endif

	std::vector<double> res(params_num, 0.0);
	x = solveLUP(B, b, params_num);
	for (int i = 0; i < params_num; i++) {
		for (int j = 0; j < params_num; j++) {
			res[i] += B[i * params_num + j] * x[j];
		}
		printf("b[i] = %.3le \t res[i] = %.3le \t error = %.3le\n", b[i], res[i], (b[i] - res[i]) / b[i]);
	}

#ifdef DEBUG 
	printf("\n   X solve\n");
	for (int i = 0; i < params_num; i++) {
		printf("%.3f ", x[i]);
	}
	printf("\n\n");
#endif
	for (int i = 0; i < params_num; i++) {
		params[i] -= x[i];
	}

	return params;
}


void read_file_into_vector(const std::string& filename, std::vector<double>& vector) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Ошибка при открытии файла: " << filename << std::endl;
		return;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::stringstream ss(line);
		double value;
		while (ss >> value) {
			vector.push_back(value);
		}
	}

	file.close();
}


int main() {
	std::vector<double> observations,
		params,
		bufers1(params_num + params_num * states_num),
		bufers2(params_num + params_num * states_num),
		bufers3(params_num + params_num * states_num),
		bufers4(params_num + params_num * states_num);
	std::vector<std::vector<double> > bufers = { bufers1, bufers2, bufers3, bufers4 };

	read_file_into_vector("Observations.txt", observations);
	// for (int i = 0; i < observations.size(); i++) {
	// 	std::cout << observations[i] << " ";
	// }
	// std::cout << "\n" << std::endl;
	read_file_into_vector("Initial parameters.txt", params);
	// for (int i = 0; i < params.size(); i++) {
	// 	std::cout << params[i] << " ";
	// }
	// std::cout << std::endl;

	for (int i = 0; i < 10; i++) {
		printf("BH mass = %.16le\n", params[18]);
		params = Gauss_Newton(params, observations, bufers);
		printf("BH mass = %.16le\n", params[18]);
		// std::cout << " " << std::endl;
	}


	/*std::vector<double> result = F(dxdP);

	for (int i = 0; i < 19; i++) {
		printf("%.3f ", result[i]);
	}
	for (int i = 0; i < 18; ++i) {
		for (int j = 0; j < 19; ++j) {
			printf("%.3f ", result[i * 19 + j + 19]);
		}
		printf("\n\n");
	}
	for (int i = 0; i < 19; i++) {
		printf("%.3f ", result[i]);
	}*/

	return 0;
}
