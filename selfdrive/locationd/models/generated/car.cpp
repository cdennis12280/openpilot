#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1229659804743780372) {
   out_1229659804743780372[0] = delta_x[0] + nom_x[0];
   out_1229659804743780372[1] = delta_x[1] + nom_x[1];
   out_1229659804743780372[2] = delta_x[2] + nom_x[2];
   out_1229659804743780372[3] = delta_x[3] + nom_x[3];
   out_1229659804743780372[4] = delta_x[4] + nom_x[4];
   out_1229659804743780372[5] = delta_x[5] + nom_x[5];
   out_1229659804743780372[6] = delta_x[6] + nom_x[6];
   out_1229659804743780372[7] = delta_x[7] + nom_x[7];
   out_1229659804743780372[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2867377941913408996) {
   out_2867377941913408996[0] = -nom_x[0] + true_x[0];
   out_2867377941913408996[1] = -nom_x[1] + true_x[1];
   out_2867377941913408996[2] = -nom_x[2] + true_x[2];
   out_2867377941913408996[3] = -nom_x[3] + true_x[3];
   out_2867377941913408996[4] = -nom_x[4] + true_x[4];
   out_2867377941913408996[5] = -nom_x[5] + true_x[5];
   out_2867377941913408996[6] = -nom_x[6] + true_x[6];
   out_2867377941913408996[7] = -nom_x[7] + true_x[7];
   out_2867377941913408996[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6562238736201300794) {
   out_6562238736201300794[0] = 1.0;
   out_6562238736201300794[1] = 0;
   out_6562238736201300794[2] = 0;
   out_6562238736201300794[3] = 0;
   out_6562238736201300794[4] = 0;
   out_6562238736201300794[5] = 0;
   out_6562238736201300794[6] = 0;
   out_6562238736201300794[7] = 0;
   out_6562238736201300794[8] = 0;
   out_6562238736201300794[9] = 0;
   out_6562238736201300794[10] = 1.0;
   out_6562238736201300794[11] = 0;
   out_6562238736201300794[12] = 0;
   out_6562238736201300794[13] = 0;
   out_6562238736201300794[14] = 0;
   out_6562238736201300794[15] = 0;
   out_6562238736201300794[16] = 0;
   out_6562238736201300794[17] = 0;
   out_6562238736201300794[18] = 0;
   out_6562238736201300794[19] = 0;
   out_6562238736201300794[20] = 1.0;
   out_6562238736201300794[21] = 0;
   out_6562238736201300794[22] = 0;
   out_6562238736201300794[23] = 0;
   out_6562238736201300794[24] = 0;
   out_6562238736201300794[25] = 0;
   out_6562238736201300794[26] = 0;
   out_6562238736201300794[27] = 0;
   out_6562238736201300794[28] = 0;
   out_6562238736201300794[29] = 0;
   out_6562238736201300794[30] = 1.0;
   out_6562238736201300794[31] = 0;
   out_6562238736201300794[32] = 0;
   out_6562238736201300794[33] = 0;
   out_6562238736201300794[34] = 0;
   out_6562238736201300794[35] = 0;
   out_6562238736201300794[36] = 0;
   out_6562238736201300794[37] = 0;
   out_6562238736201300794[38] = 0;
   out_6562238736201300794[39] = 0;
   out_6562238736201300794[40] = 1.0;
   out_6562238736201300794[41] = 0;
   out_6562238736201300794[42] = 0;
   out_6562238736201300794[43] = 0;
   out_6562238736201300794[44] = 0;
   out_6562238736201300794[45] = 0;
   out_6562238736201300794[46] = 0;
   out_6562238736201300794[47] = 0;
   out_6562238736201300794[48] = 0;
   out_6562238736201300794[49] = 0;
   out_6562238736201300794[50] = 1.0;
   out_6562238736201300794[51] = 0;
   out_6562238736201300794[52] = 0;
   out_6562238736201300794[53] = 0;
   out_6562238736201300794[54] = 0;
   out_6562238736201300794[55] = 0;
   out_6562238736201300794[56] = 0;
   out_6562238736201300794[57] = 0;
   out_6562238736201300794[58] = 0;
   out_6562238736201300794[59] = 0;
   out_6562238736201300794[60] = 1.0;
   out_6562238736201300794[61] = 0;
   out_6562238736201300794[62] = 0;
   out_6562238736201300794[63] = 0;
   out_6562238736201300794[64] = 0;
   out_6562238736201300794[65] = 0;
   out_6562238736201300794[66] = 0;
   out_6562238736201300794[67] = 0;
   out_6562238736201300794[68] = 0;
   out_6562238736201300794[69] = 0;
   out_6562238736201300794[70] = 1.0;
   out_6562238736201300794[71] = 0;
   out_6562238736201300794[72] = 0;
   out_6562238736201300794[73] = 0;
   out_6562238736201300794[74] = 0;
   out_6562238736201300794[75] = 0;
   out_6562238736201300794[76] = 0;
   out_6562238736201300794[77] = 0;
   out_6562238736201300794[78] = 0;
   out_6562238736201300794[79] = 0;
   out_6562238736201300794[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8103601707352850107) {
   out_8103601707352850107[0] = state[0];
   out_8103601707352850107[1] = state[1];
   out_8103601707352850107[2] = state[2];
   out_8103601707352850107[3] = state[3];
   out_8103601707352850107[4] = state[4];
   out_8103601707352850107[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8103601707352850107[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8103601707352850107[7] = state[7];
   out_8103601707352850107[8] = state[8];
}
void F_fun(double *state, double dt, double *out_7112769282958600023) {
   out_7112769282958600023[0] = 1;
   out_7112769282958600023[1] = 0;
   out_7112769282958600023[2] = 0;
   out_7112769282958600023[3] = 0;
   out_7112769282958600023[4] = 0;
   out_7112769282958600023[5] = 0;
   out_7112769282958600023[6] = 0;
   out_7112769282958600023[7] = 0;
   out_7112769282958600023[8] = 0;
   out_7112769282958600023[9] = 0;
   out_7112769282958600023[10] = 1;
   out_7112769282958600023[11] = 0;
   out_7112769282958600023[12] = 0;
   out_7112769282958600023[13] = 0;
   out_7112769282958600023[14] = 0;
   out_7112769282958600023[15] = 0;
   out_7112769282958600023[16] = 0;
   out_7112769282958600023[17] = 0;
   out_7112769282958600023[18] = 0;
   out_7112769282958600023[19] = 0;
   out_7112769282958600023[20] = 1;
   out_7112769282958600023[21] = 0;
   out_7112769282958600023[22] = 0;
   out_7112769282958600023[23] = 0;
   out_7112769282958600023[24] = 0;
   out_7112769282958600023[25] = 0;
   out_7112769282958600023[26] = 0;
   out_7112769282958600023[27] = 0;
   out_7112769282958600023[28] = 0;
   out_7112769282958600023[29] = 0;
   out_7112769282958600023[30] = 1;
   out_7112769282958600023[31] = 0;
   out_7112769282958600023[32] = 0;
   out_7112769282958600023[33] = 0;
   out_7112769282958600023[34] = 0;
   out_7112769282958600023[35] = 0;
   out_7112769282958600023[36] = 0;
   out_7112769282958600023[37] = 0;
   out_7112769282958600023[38] = 0;
   out_7112769282958600023[39] = 0;
   out_7112769282958600023[40] = 1;
   out_7112769282958600023[41] = 0;
   out_7112769282958600023[42] = 0;
   out_7112769282958600023[43] = 0;
   out_7112769282958600023[44] = 0;
   out_7112769282958600023[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7112769282958600023[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7112769282958600023[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7112769282958600023[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7112769282958600023[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7112769282958600023[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7112769282958600023[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7112769282958600023[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7112769282958600023[53] = -9.8000000000000007*dt;
   out_7112769282958600023[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7112769282958600023[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7112769282958600023[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7112769282958600023[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7112769282958600023[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7112769282958600023[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7112769282958600023[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7112769282958600023[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7112769282958600023[62] = 0;
   out_7112769282958600023[63] = 0;
   out_7112769282958600023[64] = 0;
   out_7112769282958600023[65] = 0;
   out_7112769282958600023[66] = 0;
   out_7112769282958600023[67] = 0;
   out_7112769282958600023[68] = 0;
   out_7112769282958600023[69] = 0;
   out_7112769282958600023[70] = 1;
   out_7112769282958600023[71] = 0;
   out_7112769282958600023[72] = 0;
   out_7112769282958600023[73] = 0;
   out_7112769282958600023[74] = 0;
   out_7112769282958600023[75] = 0;
   out_7112769282958600023[76] = 0;
   out_7112769282958600023[77] = 0;
   out_7112769282958600023[78] = 0;
   out_7112769282958600023[79] = 0;
   out_7112769282958600023[80] = 1;
}
void h_25(double *state, double *unused, double *out_8606043731816012267) {
   out_8606043731816012267[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1177661407255701104) {
   out_1177661407255701104[0] = 0;
   out_1177661407255701104[1] = 0;
   out_1177661407255701104[2] = 0;
   out_1177661407255701104[3] = 0;
   out_1177661407255701104[4] = 0;
   out_1177661407255701104[5] = 0;
   out_1177661407255701104[6] = 1;
   out_1177661407255701104[7] = 0;
   out_1177661407255701104[8] = 0;
}
void h_24(double *state, double *unused, double *out_5898819030352362688) {
   out_5898819030352362688[0] = state[4];
   out_5898819030352362688[1] = state[5];
}
void H_24(double *state, double *unused, double *out_999553016351448869) {
   out_999553016351448869[0] = 0;
   out_999553016351448869[1] = 0;
   out_999553016351448869[2] = 0;
   out_999553016351448869[3] = 0;
   out_999553016351448869[4] = 1;
   out_999553016351448869[5] = 0;
   out_999553016351448869[6] = 0;
   out_999553016351448869[7] = 0;
   out_999553016351448869[8] = 0;
   out_999553016351448869[9] = 0;
   out_999553016351448869[10] = 0;
   out_999553016351448869[11] = 0;
   out_999553016351448869[12] = 0;
   out_999553016351448869[13] = 0;
   out_999553016351448869[14] = 1;
   out_999553016351448869[15] = 0;
   out_999553016351448869[16] = 0;
   out_999553016351448869[17] = 0;
}
void h_30(double *state, double *unused, double *out_4744360255387510615) {
   out_4744360255387510615[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1340671551251547523) {
   out_1340671551251547523[0] = 0;
   out_1340671551251547523[1] = 0;
   out_1340671551251547523[2] = 0;
   out_1340671551251547523[3] = 0;
   out_1340671551251547523[4] = 1;
   out_1340671551251547523[5] = 0;
   out_1340671551251547523[6] = 0;
   out_1340671551251547523[7] = 0;
   out_1340671551251547523[8] = 0;
}
void h_26(double *state, double *unused, double *out_740723775953321840) {
   out_740723775953321840[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4919164726129757328) {
   out_4919164726129757328[0] = 0;
   out_4919164726129757328[1] = 0;
   out_4919164726129757328[2] = 0;
   out_4919164726129757328[3] = 0;
   out_4919164726129757328[4] = 0;
   out_4919164726129757328[5] = 0;
   out_4919164726129757328[6] = 0;
   out_4919164726129757328[7] = 1;
   out_4919164726129757328[8] = 0;
}
void h_27(double *state, double *unused, double *out_7623959638480714155) {
   out_7623959638480714155[0] = state[3];
}
void H_27(double *state, double *unused, double *out_834091760548877388) {
   out_834091760548877388[0] = 0;
   out_834091760548877388[1] = 0;
   out_834091760548877388[2] = 0;
   out_834091760548877388[3] = 1;
   out_834091760548877388[4] = 0;
   out_834091760548877388[5] = 0;
   out_834091760548877388[6] = 0;
   out_834091760548877388[7] = 0;
   out_834091760548877388[8] = 0;
}
void h_29(double *state, double *unused, double *out_4051949208878486525) {
   out_4051949208878486525[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1850902895565939707) {
   out_1850902895565939707[0] = 0;
   out_1850902895565939707[1] = 1;
   out_1850902895565939707[2] = 0;
   out_1850902895565939707[3] = 0;
   out_1850902895565939707[4] = 0;
   out_1850902895565939707[5] = 0;
   out_1850902895565939707[6] = 0;
   out_1850902895565939707[7] = 0;
   out_1850902895565939707[8] = 0;
}
void h_28(double *state, double *unused, double *out_8791703439556855892) {
   out_8791703439556855892[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3231496121503590867) {
   out_3231496121503590867[0] = 1;
   out_3231496121503590867[1] = 0;
   out_3231496121503590867[2] = 0;
   out_3231496121503590867[3] = 0;
   out_3231496121503590867[4] = 0;
   out_3231496121503590867[5] = 0;
   out_3231496121503590867[6] = 0;
   out_3231496121503590867[7] = 0;
   out_3231496121503590867[8] = 0;
}
void h_31(double *state, double *unused, double *out_4862367649320887359) {
   out_4862367649320887359[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5545372828363108804) {
   out_5545372828363108804[0] = 0;
   out_5545372828363108804[1] = 0;
   out_5545372828363108804[2] = 0;
   out_5545372828363108804[3] = 0;
   out_5545372828363108804[4] = 0;
   out_5545372828363108804[5] = 0;
   out_5545372828363108804[6] = 0;
   out_5545372828363108804[7] = 0;
   out_5545372828363108804[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1229659804743780372) {
  err_fun(nom_x, delta_x, out_1229659804743780372);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2867377941913408996) {
  inv_err_fun(nom_x, true_x, out_2867377941913408996);
}
void car_H_mod_fun(double *state, double *out_6562238736201300794) {
  H_mod_fun(state, out_6562238736201300794);
}
void car_f_fun(double *state, double dt, double *out_8103601707352850107) {
  f_fun(state,  dt, out_8103601707352850107);
}
void car_F_fun(double *state, double dt, double *out_7112769282958600023) {
  F_fun(state,  dt, out_7112769282958600023);
}
void car_h_25(double *state, double *unused, double *out_8606043731816012267) {
  h_25(state, unused, out_8606043731816012267);
}
void car_H_25(double *state, double *unused, double *out_1177661407255701104) {
  H_25(state, unused, out_1177661407255701104);
}
void car_h_24(double *state, double *unused, double *out_5898819030352362688) {
  h_24(state, unused, out_5898819030352362688);
}
void car_H_24(double *state, double *unused, double *out_999553016351448869) {
  H_24(state, unused, out_999553016351448869);
}
void car_h_30(double *state, double *unused, double *out_4744360255387510615) {
  h_30(state, unused, out_4744360255387510615);
}
void car_H_30(double *state, double *unused, double *out_1340671551251547523) {
  H_30(state, unused, out_1340671551251547523);
}
void car_h_26(double *state, double *unused, double *out_740723775953321840) {
  h_26(state, unused, out_740723775953321840);
}
void car_H_26(double *state, double *unused, double *out_4919164726129757328) {
  H_26(state, unused, out_4919164726129757328);
}
void car_h_27(double *state, double *unused, double *out_7623959638480714155) {
  h_27(state, unused, out_7623959638480714155);
}
void car_H_27(double *state, double *unused, double *out_834091760548877388) {
  H_27(state, unused, out_834091760548877388);
}
void car_h_29(double *state, double *unused, double *out_4051949208878486525) {
  h_29(state, unused, out_4051949208878486525);
}
void car_H_29(double *state, double *unused, double *out_1850902895565939707) {
  H_29(state, unused, out_1850902895565939707);
}
void car_h_28(double *state, double *unused, double *out_8791703439556855892) {
  h_28(state, unused, out_8791703439556855892);
}
void car_H_28(double *state, double *unused, double *out_3231496121503590867) {
  H_28(state, unused, out_3231496121503590867);
}
void car_h_31(double *state, double *unused, double *out_4862367649320887359) {
  h_31(state, unused, out_4862367649320887359);
}
void car_H_31(double *state, double *unused, double *out_5545372828363108804) {
  H_31(state, unused, out_5545372828363108804);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
