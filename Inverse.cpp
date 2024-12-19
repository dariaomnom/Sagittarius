#include <iostream>
#include <cmath>
#include <vector>

const int stars_num = 3;
const int states_num = 6 * 3; // (3 coords + 3 velocities) * 3 stars
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

double RA(std::vector<double>& params) {
	return 0;
}

double Decl(std::vector<double>& params) {
	return 0;
}

void RK4(std::vector<double>& params, std::vector<double>& stateBuf, std::vector<double>& k1, std::vector<double>& k2, std::vector<double>& k3, double h)
{
	k1 = F_GR_3(params);
	for (int i = 0; i < states_num; i++) {
		stateBuf[i] = params[i] + 0.5 * h * k1[i];
	}

	k2 = F_GR_3(params);
	for (int i = 0; i < states_num; i++) {
		stateBuf[i] = params[i] + 0.75 * h * k2[i];
	}

	k3 = F_GR_3(params);
	for (int i = 0; i < states_num; i++) {
		params[i] += h * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]);
	}
}

std::vector<double> residuals(std::vector<double>& params, std::vector<double>& observations, std::vector<std::vector<double>>& bufers) {
	std::vector<double> residuals(observation_num);
	double t = start_time;
	int observ_counter = 0;
	while (t < end_time) {
		for (int i = 0; i < step / h; i++, t += h) {
			RK4(params, bufers[0], bufers[1], bufers[2], bufers[3], h);
		}
		if (t == observations[observ_counter * 3]) {
			residuals[observ_counter * 2] = observations[observ_counter * 3 + 1] - RA(params);
			residuals[observ_counter * 2 + 1] = observations[observ_counter * 3 + 2] - Decl(params);
		}
	}



	return residuals;
}

int main() {

	return 0;
}