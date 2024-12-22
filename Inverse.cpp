#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

const int stars_num = 3; // S2, S38, S55
const int states_num = 6 * 3; // (3 coords + 3 velocities) * 3 stars
const int params_num = states_num + 1; // states + mass BH
const int observation_num = 2 * (33 + 29 + 25); // alpha, delta for each observation of three stars
const double RAD_TO_ARCSECONDS = 206264.8;
const double C = 299792458.0;
const double start_time = 2002.500;
const double end_time = 2016;
const double Y = 31558149.0;
const double STEP = Y / 10000.0;
const std::vector<double> BH_cart = std::vector<double>();
const std::vector<double> BH_eq = std::vector<double>();
std::vector<std::vector<double> > observations;

std::vector<double> F_GR_3(std::vector<double>& condition) {
	std::vector<double> d_condition(states_num), R(stars_num), Vsqrd(stars_num), VxR(stars_num);
	double GM = condition.back();

	R[0] = sqrt(condition[0] * condition[0] + condition[1] * condition[1] + condition[2] * condition[2]);
	R[1] = sqrt(condition[3] * condition[3] + condition[4] * condition[4] + condition[5] * condition[5]);
	R[2] = sqrt(condition[6] * condition[6] + condition[7] * condition[7] + condition[8] * condition[8]);

	Vsqrd[0] = condition[9] * condition[9] + condition[10] * condition[10] + condition[11] * condition[11];
	Vsqrd[1] = condition[12] * condition[12] + condition[13] * condition[13] + condition[14] * condition[14];
	Vsqrd[2] = condition[15] * condition[15] + condition[16] * condition[16] + condition[17] * condition[17];

	VxR[0] = condition[0] * condition[9] + condition[1] * condition[10] + condition[2] * condition[11];
	VxR[1] = condition[3] * condition[12] + condition[4] * condition[13] + condition[5] * condition[14];
	VxR[2] = condition[6] * condition[15] + condition[7] * condition[16] + condition[8] * condition[17];

	for (int i = 0; i < states_num / 2; i++) {
		d_condition[i] = condition[i + states_num / 2];
	}

	for (int i = 0; i < states_num / 2; i++) {
		d_condition[i + states_num / 2] = -GM / (C * C * R[i / 3] * R[i / 3] * R[i / 3]) \
			* (condition[i] * (C * C - 4 * GM / R[i / 3] + Vsqrd[i / 3]) - 4 * condition[i + states_num / 2] * VxR[i / 3]);
	}


	std::vector<double> result(states_num * params_num, 0.0), dFdx(states_num * states_num, 0.0), dFdP(states_num / 2, 0.0);



	for (int i = 0; i < states_num / 2; i++) {
		dFdx[states_num * i + i + states_num / 2] = 1;
	}

	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				int ik = 9 + i + 3 * k;
				int jk = j + 3 * k;
				if (ik - 9 == jk)
					dFdx[states_num * ik + jk] = -M * (pow(R[k], 2) - 3 * pow(condition[jk], 2)) / pow(R[k], 5);
				else
					dFdx[states_num * ik + jk] = M * 3 * condition[ik - 9] * condition[jk] / pow(R[k], 5);
			}
		}
	}

	// result sum of martices
	for (int i = 0; i < states_num; i++) { // (dF/dP) + (dF/dx)(dx/dP)
		for (int j = 0; j < params_num; j++) {

			if (j == params_num - 1 && i >= states_num / 2)
				result[i * params_num + j] = dFdP[i - states_num / 2];

			for (int k = 0; k < states_num; k++) {
				result[i * params_num + j] += dFdx[i * states_num + k] * dxdP[k * params_num + j];
			}

		}
	}
	return result;
}


double gRA(std::vector<double>& params, int observ_counter) {
	int k;
	if (observ_counter <= 33) k = 0;
	else if (observ_counter <= 62) k = 1;
	else k = 2;

	double x = params[3 * k] + BH_cart[0];
	double y = params[3 * k + 1] + BH_cart[1];
	double z = params[3 * k + 2] + BH_cart[2];

	double alpha = atan2(y, x);

	return alpha * RAD_TO_ARCSECONDS - BH_eq[0];
}


double gDecl(std::vector<double>& params, int observ_counter) {
	int k;
	if (observ_counter <= 33) k = 0;
	else if (observ_counter <= 62) k = 1;
	else k = 2;

	double x = params[3 * k] + BH_cart[0];
	double y = params[3 * k + 1] + BH_cart[1];
	double z = params[3 * k + 2] + BH_cart[2];

	double r = sqrt(x * x + y * y + z * z);
	double delta = asin(z / r);

	return delta * RAD_TO_ARCSECONDS - BH_eq[1];
}


std::vector<double> dRA(std::vector<double>& params) {
	std::vector<double> result;
	return result;
}


std::vector<double> dDecl(std::vector<double>& params) {
	std::vector<double> result;
	return result;
}


void RK4_step(std::vector<double>& condition, std::vector<std::vector<double> >& bufers, double h)
{
	std::vector<double> stateBuf = bufers[0], k1 = bufers[1], k2 = bufers[2], k3 = bufers[3];
	k1 = F_GR_3(condition);
	for (int i = 0; i < states_num; i++) {
		stateBuf[i] = condition[i] + 0.5 * h * k1[i];
	}

	k2 = F_GR_3(state);
	for (int i = 0; i < states_num; i++) {
		stateBuf[i] = state[i] + 0.75 * h * k2[i];
	}

	k3 = F_GR_3(state);
	for (int i = 0; i < states_num; i++) {
		state[i] += h * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]);
	}
}


std::vector<double> Gauss_Newton(std::vector<double>& params, std::vector<double>& observations, std::vector<std::vector<double> >& bufers) {
	std::vector<double> dgdx(states_num, 0.0), residuals(observation_num, 0.0), params_new(params_num, 0.0), condition(params_num + states_num * params_num);
	for (int i = 0; i < params_num; i++) condition[i] = params[i];
	for (int i = 0; i < states_num * params_num; i++) {
		if (i / params_num == i % params_num) condition[i + params_num] = 1.0;
	}
	double t = start_time;
	double step = STEP;
	int observ_counter = 0;
	while (t < end_time) {
		RK4_step(condition, bufers, step);

		if (observations[observ_counter * 5] - t < step) {
			step = observations[observ_counter * 5] - t;
			RK4_step(condition, bufers, step);
			residuals[observ_counter * 2] = observations[observ_counter * 5 + 1] - gRA(params, observ_counter);
			residuals[observ_counter * 2 + 1] = observations[observ_counter * 5 + 2] - gDecl(params, observ_counter);

			// Заполнение соответствующих строк матрицы А

			observ_counter++;
			step = STEP;
		}
	}
	// params_new = params - (AT * W * A)^-1 * AT * W * residuals
	return params_new;
}


void read_file_into_matrix(const std::string& filename, std::vector<std::vector<double> >& matrix) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Ошибка при открытии файла: " << filename << std::endl;
		return;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::vector<double> row;
		std::stringstream ss(line);
		double value;

		while (ss >> value) {
			row.push_back(value);
		}

		matrix.push_back(row);
	}

	file.close();
}

void read_observations(std::vector<std::vector<double> >& observations) {
	observations.clear();
	read_file_into_matrix("S2_Article.txt", observations);
	read_file_into_matrix("S38_Article.txt", observations);
	read_file_into_matrix("S55_Article.txt", observations);
}

int main() {
	read_observations(observations);
	// test read
	// for (size_t i = 0; i < observations.size(); ++i) {
	//     for (size_t j = 0; j < observations[i].size(); ++j) {
	//         std::cout << observations[i][j] << " ";
	//     }
	//     std::cout << std::endl;
	// }
	std::vector<double> dxdP(18 * 19, 0.0); // Матрица 18x19, заполнена нулями
	for (int i = 0; i < 18; ++i) {
		dxdP[i * 19 + i] = 1.0; // Установка единиц на главной диагонали (18x18 часть)
	}

	std::vector<double> result = F_dxdP(dxdP, state);

	// std::cout << "Матрица dxdP (18x19):" << std::endl;
	// for (int i = 0; i < 18; ++i) {
	//     for (int j = 0; j < 19; ++j) {
	//         std::cout << std::setw(2) << dxdP[i * 19 + j] << " ";
	//     }
	//     std::cout << std::endl;
	// }

	return 0;
}
