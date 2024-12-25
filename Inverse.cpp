#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#define DEBUG 1
const int stars_num = 3; // S2, S38, S55
const int states_num = 6 * 3; // (3 coords + 3 velocities) * 3 stars
const int params_num = states_num + 1; // states + mass BH
const int observation_num = 33 + 29 + 26; // alpha, delta for each observation of three stars
const double RAD_TO_ARCSECONDS = 206264.8;
const double C = 299792458.0;
const double start_time = 2002.0;
const double end_time = 2016.0;
const double Y = 31558149.0;
const double STEP = Y / 10000.0;
const std::vector<double> BH_cart = { -1.4388065248297531e+19, -2.2977127354968464e+20, -1.2758089026580264e+20 };
const std::vector<double> BH_eq = { 959100.7014, -104377.4682 };

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

	k2 = F(condition);
	for (int i = 0; i < params_num + states_num * params_num; i++) {
		conditionBuf[i] = condition[i] + 0.75 * h * k2[i];
	}

	k3 = F(condition);
	for (int i = 0; i < params_num + states_num * params_num; i++) {
		condition[i] += h * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]);
	}
}


double g(std::vector<double>& condition, int observ_counter, int k) { // testing

	double x = condition[3 * k] + BH_cart[0];
	double y = condition[3 * k + 1] + BH_cart[1];
	double z = condition[3 * k + 2] + BH_cart[2];
	double r = sqrt(x * x + y * y + z * z);

	if (observ_counter / observation_num == 0) {
		double alpha = atan2(y, x);
		return alpha * RAD_TO_ARCSECONDS - BH_eq[0];
	}
	if (observ_counter / observation_num == 1) {
		double delta = asin(z / r);
		return delta * RAD_TO_ARCSECONDS - BH_eq[1];
	}
	return 0;
}


std::vector<double> dgdx(std::vector<double>& condition, int observ_counter, int k) { // testing
	std::vector<double> result(states_num, 0.0);

	double x = condition[3 * k] + BH_cart[0];
	double y = condition[3 * k + 1] + BH_cart[1];
	double z = condition[3 * k + 2] + BH_cart[2];
	double r = sqrt(x * x + y * y + z * z);

	if (observ_counter / observation_num == 0) {
		result[k] = -y / (x * x + y * y) * RAD_TO_ARCSECONDS;
		result[k + 1] = x / (x * x + y * y) * RAD_TO_ARCSECONDS;
	}
	if (observ_counter / observation_num == 1) {
		result[k] = -z * x / (pow(r, 3.0) * pow(1 - pow(z / r, 2.0), 1.0 / 2)) * RAD_TO_ARCSECONDS;
		result[k + 1] = -z * y / (pow(r, 3.0) * pow(1 - pow(z / r, 2.0), 1.0 / 2)) * RAD_TO_ARCSECONDS;
		result[k + 2] = (x * x + y * y) / (pow(r, 3.0) * pow(1 - pow(z / r, 2.0), 1.0 / 2)) * RAD_TO_ARCSECONDS;
	}

	return result;
}


std::vector<double> CholeskyDecomposition(const std::vector<double>& A, int n) {
	std::vector<double> L(n * n, 0.0);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= i; ++j) {
			double sum = 0.0;

			if (j == i) {
				for (int k = 0; k < j; ++k)
					sum += L[j * n + k] * L[j * n + k];
				L[j * n + j] = std::sqrt(A[j * n + j] - sum);
			}
			else {
				for (int k = 0; k < j; ++k)
					sum += L[i * n + k] * L[j * n + k];
				L[i * n + j] = (A[i * n + j] - sum) / L[j * n + j];
			}
		}
	}

	return L;
}

// Прямой ход для L * y = B
std::vector<double> ForwardSubstitution(const std::vector<double>& L, const std::vector<double>& B, int n) {
	std::vector<double> y(n, 0.0);

	for (int i = 0; i < n; ++i) {
		double sum = 0.0;
		for (int j = 0; j < i; ++j)
			sum += L[i * n + j] * y[j];
		y[i] = (B[i] - sum) / L[i * n + i];
	}

	return y;
}

// Обратный ход для L^T * x = y
std::vector<double> BackwardSubstitution(const std::vector<double>& L, const std::vector<double>& y, int n) {
	std::vector<double> x(n, 0.0);

	for (int i = n - 1; i >= 0; --i) {
		double sum = 0.0;
		for (int j = i + 1; j < n; ++j)
			sum += L[j * n + i] * x[j];
		x[i] = (y[i] - sum) / L[i * n + i];
	}

	return x;
}

// Решение системы AX = B с использованием разложения Холецкого
std::vector<double> SolveCholesky(const std::vector<double>& A, const std::vector<double>& B, int n) {

	// Разложение A = L * L^T
	std::vector<double> L = CholeskyDecomposition(A, n);

	// Решение L * y = B
	std::vector<double> y = ForwardSubstitution(L, B, n);

	// Решение L^T * x = y
	std::vector<double> x = BackwardSubstitution(L, y, n);

	return x;
}

std::vector<double> Gauss_Newton(std::vector<double>& params, std::vector<double>& observations, std::vector<std::vector<double> >& bufers) {
	std::vector<double> condition(params_num + states_num * params_num),
		b(params_num, 0.0),                                 // AT * Wsqrd * r(beta)
		x(params_num, 0.0),                                 // AT * Wsqrd * r(beta)
		B(params_num * params_num, 0.0),                    // (AT * W * A)
		dg(states_num, 0.0),
		Ak(params_num, 0.0);

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
			RK4_step(condition, bufers, step);
			t += step;
			if (DEBUG) {
				std::cout << "\n";
				std::cout << observ_counter << std::endl;
				std::cout << observations[observ_counter * 6] << std::endl;
				std::cout << "step " << step << std::endl;
			}


			// RA
			dg = dgdx(condition, observ_counter, observations[observ_counter * 6 + 5]);
			for (int i = 0; i < params_num; i++) {	// line of A matrix -(dg/dx * dx/dP)_k
				Ak[i] = 0.0;
				for (int j = 0; j < states_num; j++) {
					Ak[i] -= dg[j] * condition[j * params_num + i + params_num];
				}
			}

			for (int i = 0; i < params_num; i++) {	// A * Wsqrd * r
				b[i] += Ak[i] / observations[observ_counter * 6 + 3] * (observations[observ_counter * 6 + 1] - g(condition, observ_counter, observations[observ_counter * 6 + 5]));
			}

			for (int i = 0; i < params_num; i++) {	// AT * W * A
				for (int j = 0; j < params_num; j++) {
					B[i * params_num + j] += Ak[i] * Ak[j] / pow(observations[observ_counter * 6 + 3], 2.0);
				}
			}

			// Decl
			dg = dgdx(condition, observ_counter, observations[observ_counter * 6 + 5]);
			for (int i = 0; i < params_num; i++) {	// line of A matrix -(dg/dx * dx/dP)_k
				Ak[i] = 0.0;
				for (int j = 0; j < states_num; j++) {
					Ak[i] -= dg[j] * condition[j * params_num + i + params_num];
				}
			}

			for (int i = 0; i < params_num; i++) {	// A * Wsqrd * r
				b[i] += Ak[i] / observations[observ_counter * 6 + 4] * (observations[observ_counter * 6 + 2] - g(condition, observ_counter, observations[observ_counter * 6 + 5]));
			}

			for (int i = 0; i < params_num; i++) {	// AT * W * A
				for (int j = 0; j < params_num; j++) {
					B[i * params_num + j] += Ak[i] * Ak[j] / pow(observations[observ_counter * 6 + 4], 2.0);
				}
			}

			observ_counter++;
			step = STEP;
			/*if (DEBUG) {
				for (int i = 0; i < params_num; i++) {
					for (int j = 0; j < params_num; j++) {
						printf("%.3f ", B[i * params_num + j]);
					}
					printf("\n");
				}
			}*/
			/*if (DEBUG) {
				for (int i = 0; i < params_num; i++) {
					printf("%.3f ", b[i]);
				}
				printf("\n");
			}*/
		}
		else {
			RK4_step(condition, bufers, step);
			t += step;
		}

	}
	if (DEBUG) {
		for (int i = 0; i < params_num; i++) {
			printf("%.3f ", b[i]);
		}
		printf("\n\n");
	}
	x = SolveCholesky(B, b, params_num);
	if (DEBUG) {
		for (int i = 0; i < params_num; i++) {
			printf("%.3f ", x[i]);
		}
		printf("\n\n");
	}
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
