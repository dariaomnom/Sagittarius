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
const double step = Y / 1000.0;
double const h = Y / 10000.0;
const std::vector<double> BH_cart = std::vector<double>();
const std::vector<double> BH_eq = std::vector<double>();
std::vector<std::vector<double> > observations;

std::vector<double> F_GR_3(std::vector<double>& state) {
	std::vector<double> d_state(states_num), R(stars_num), Vsqrd(stars_num), VxR(stars_num);
	double GM = state.back();

	for (int i = 0; i < states_num / 2; i++) {
		d_state[i] = state[i + states_num / 2];
	}

	R[0] = sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
	R[1] = sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
	R[2] = sqrt(state[6] * state[6] + state[7] * state[7] + state[8] * state[8]);

	Vsqrd[0] = state[9] * state[9] + state[10] * state[10] + state[11] * state[11];
	Vsqrd[1] = state[12] * state[12] + state[13] * state[13] + state[14] * state[14];
	Vsqrd[2] = state[15] * state[15] + state[16] * state[16] + state[17] * state[17];

	VxR[0] = state[0] * state[9] + state[1] * state[10] + state[2] * state[11];
	VxR[1] = state[3] * state[12] + state[4] * state[13] + state[5] * state[14];
	VxR[2] = state[6] * state[15] + state[7] * state[16] + state[8] * state[17];

	for (int i = 0; i < states_num / 2; i++) {
		d_state[i + states_num / 2] = -GM / (C * C * R[i / 3] * R[i / 3] * R[i / 3]) * (state[i] * (C * C - 4 * GM / R[i / 3] + Vsqrd[i / 3]) - 4 * state[i + states_num / 2] * VxR[i / 3]);
	}
	return d_state;
}

std::vector<double> F_dxdP(std::vector<double>& dxdP, std::vector<double>& state) { // testing
	std::vector<double> result(states_num * params_num, 0.0), dFdx(states_num * states_num, 0.0), dFdP(states_num / 2, 0.0);
	double M = state.back();
	std::vector<double> R(stars_num);
	R[0] = sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
	R[1] = sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
	R[2] = sqrt(state[6] * state[6] + state[7] * state[7] + state[8] * state[8]);

	// last column

	for (int i = 0; i < states_num / 2; i++) {
		dFdP[i] = -M * state[i] / pow(R[i / 3], 3.);
	}
	std::cout << std::fixed << std::setprecision(3);
    std::cout << "Последний столбец:" << std::endl;
    for (int i = 0; i < dFdP.size(); i++) {
		if (dFdP[i] == 0.0) {
			std::cout << '0';
		} else {
			std::cout << dFdP[i];
		}
        std::cout << std::endl;    
    }

	// E-matrix
	for (int i = 0; i < states_num / 2; i++) {
		dFdx[states_num * i + i + states_num / 2] = 1;
	}
	std::cout << "Единичная матрица:" << std::endl;
	for (int i = 0; i < states_num; i++) {
		for (int j = 0; j < states_num; j++) {
			if (dFdx[i * states_num + j] == 0.0) {
				std::cout << std::setw(7) << '0';
			} else {
				std::cout << std::setw(7) << dFdx[i * states_num + j];
			}
		}
		std::cout << std::endl;
	}
	
	// left down block
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				int ik = 9 + i + 3 * k;
				int jk = j + 3 * k;
				if (ik - 9 == jk) 
					dFdx[states_num * ik + jk] = -M * (pow(R[k], 2) - 3 * pow(state[jk], 2)) / pow(R[k], 5);
				else 
					dFdx[states_num * ik + jk] = M * 3 * state[ik - 9] * state[jk] / pow(R[k], 5);
			}
		}
	}
	std::cout << std::fixed << std::setprecision(3);
    std::cout << "Левый нижний блок:" << std::endl;
	for (int i = 0; i < states_num; i++) {
		for (int j = 0; j < states_num; j++) {
			if (dFdx[i * states_num + j] == 0.0) {
				std::cout << std::setw(7) << '0';
			} else {
				std::cout << std::setw(7) << dFdx[i * states_num + j];
			}
		}
		std::cout << std::endl;
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

void RK4_orbit(std::vector<double>& state, std::vector<std::vector<double> >& bufers, double h)
{
	std::vector<double> stateBuf = bufers[0], k1 = bufers[1], k2 = bufers[2], k3 = bufers[3];
	k1 = F_GR_3(state);
	for (int i = 0; i < states_num; i++) {
		stateBuf[i] = state[i] + 0.5 * h * k1[i];
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

void RK4_isohronic(std::vector<double>& dxdP, std::vector<double>& params, std::vector<std::vector<double> >& bufers, double h)
{
	std::vector<double> dxdPBuf = bufers[0], k1 = bufers[1], k2 = bufers[2], k3 = bufers[3];
	k1 = F_dxdP(dxdP, params);
	for (int i = 0; i < states_num * params_num; i++) {
		dxdPBuf[i] = dxdP[i] + 0.5 * h * k1[i];
	}

	k2 = F_dxdP(dxdP, params);
	for (int i = 0; i < states_num * params_num; i++) {
		dxdPBuf[i] = dxdP[i] + 0.75 * h * k2[i];
	}

	k3 = F_dxdP(dxdP, params);
	for (int i = 0; i < states_num * params_num; i++) {
		dxdP[i] += h * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]);
	}
}

std::vector<double> Gauss_Newton(std::vector<double>& params, std::vector<double>& dxdP, std::vector<double>& W, std::vector<double>& observations, std::vector<std::vector<double> >& bufers) {
	std::vector<double> dgdx(states_num), A(observation_num * params_num), residuals(observation_num), params_new(params_num);
	params_new = params;
	double t = start_time;
	int observ_counter = 0;
	while (t < end_time) {
		RK4_isohronic(dxdP, params, bufers, step);
		for (int i = 0; i < step / h; i++, t += h) {
			RK4_orbit(params, bufers, h);
		}

		if (t == observations[observ_counter * 3]) {
			residuals[observ_counter * 2] = observations[observ_counter * 3 + 1] - gRA(params, observ_counter);
			residuals[observ_counter * 2 + 1] = observations[observ_counter * 3 + 2] - gDecl(params, observ_counter);

			// Заполнение соответствующих строк матрицы А

			observ_counter++;
		}
	}
	// params_new -= (AT * W * A)^-1 * AT * W * residuals
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


	// test F_dxdP
    std::vector<double> state(params_num);
    for (int i = 0; i < params_num; ++i) {
        state[i] = static_cast<double>(i + 1);
    }

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

	std::cout << std::fixed << std::setprecision(3);
    std::cout << "Результат функции F_dxdP:" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
		if (result[i] == 0.0) {
			std::cout << std::setw(7) << '0';
		} else {
			std::cout << std::setw(7) << result[i];
		}
        if ((i + 1) % params_num == 0) {
            std::cout << std::endl;
        }
    }

	return 0;
}
