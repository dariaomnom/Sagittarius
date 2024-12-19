#include <iostream>
#include <cmath>
#include <vector>

const int stars_num = 3; // S2, S38, S55
const int states_num = 6 * 3; // (3 coords + 3 velocities) * 3 stars
const int params_num = states_num + 1; // states + mass BH
const int observation_num = 2 * (33 + 29 + 25); // alpha, delta for each observation of three stars
const double C = 299792458.0;
const double start_time = 2002.500;
const double end_time = 2016;
const double Y = 31558149.0;
const double step = Y / 1000.0;
double const h = Y / 10000.0;

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
	std::vector<double> result(states_num * params_num, 0.0), dFdP(states_num * params_num, 0.0), dFdx(states_num * states_num, 0.0);
	double M = state.back();
	std::vector<double> R(stars_num);
	R[0] = sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
	R[1] = sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
	R[2] = sqrt(state[6] * state[6] + state[7] * state[7] + state[8] * state[8]);

	for (int i = states_num / 2; i < states_num; i++) {
		dFdP[states_num * i + params_num - 1] = -M * state[i - states_num / 2] / pow(R[i / 3 - 3], 3.);
	}
	for (int i = states_num / 2; i < states_num; i++) {
		dFdx[states_num * i + i] = 1;
	}
	for (int k = 0; k < 3; k++) {
		for (int i = 9; i < 12; i++) {
			for (int j = 0; j < 3; j++) {
				int ik = i + 3 * k;
				int jk = j + 3 * k;
				if (i - 9 == j) dFdx[states_num * ik + jk] = -M * (pow(R[k], 2) - 3 * pow(state[j], 2)) / pow(R[k], 5);
				else dFdx[states_num * ik + jk] = M * 3 * state[i - 9] * state[j] / pow(R[k], 5);
			}
		}
	}
	for (int i = 0; i < states_num; i++) { // (dF/dP) + (dF/dx)(dx/dP)
		for (int j = 0; j < params_num; j++) {
			result[i * states_num + j] = dFdP[i * states_num + j];
			for (int k = 0; k < states_num; k++) {
				result[i * states_num + j] += dFdx[i * states_num + k] * dxdP[k * states_num + j];
			}
		}
	}
	return result;
}

double RA(std::vector<double>& params) {
	return 0;
}

double Decl(std::vector<double>& params) {
	return 0;
}

std::vector<double> dRA(std::vector<double>& params) {
	std::vector<double> result;
	return result;
}

std::vector<double> dDecl(std::vector<double>& params) {
	std::vector<double> result;
	return result;
}

void RK4_orbit(std::vector<double>& state, std::vector<std::vector<double>>& bufers, double h)
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

void RK4_isohronic(std::vector<double>& dxdP, std::vector<double>& params, std::vector<std::vector<double>>& bufers, double h)
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

std::vector<double> Gauss_Newton(std::vector<double>& params, std::vector<double>& dxdP, std::vector<double>& W, std::vector<double>& observations, std::vector<std::vector<double>>& bufers) {
	std::vector<double> dgdx(states_num), A(observation_num * params_num), residuals(observation_num), params_new(params_num);
	double t = start_time;
	int observ_counter = 0;
	while (t < end_time) {
		RK4_isohronic(dxdP, params, bufers, step);
		for (int i = 0; i < step / h; i++, t += h) {
			RK4_orbit(params, bufers, h);
		}

		if (t == observations[observ_counter * 3]) {
			residuals[observ_counter * 2] = observations[observ_counter * 3 + 1] - RA(params);
			residuals[observ_counter * 2 + 1] = observations[observ_counter * 3 + 2] - Decl(params);

			// Заполнение соответствующих строк матрицы А

			observ_counter++;
		}
	}
	// params_new = params - (AT * W * A)^-1 * AT * W * residuals
	return params_new;
}

int main() {

	return 0;
}