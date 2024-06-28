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
void err_fun(double *nom_x, double *delta_x, double *out_2012298189817593187) {
   out_2012298189817593187[0] = delta_x[0] + nom_x[0];
   out_2012298189817593187[1] = delta_x[1] + nom_x[1];
   out_2012298189817593187[2] = delta_x[2] + nom_x[2];
   out_2012298189817593187[3] = delta_x[3] + nom_x[3];
   out_2012298189817593187[4] = delta_x[4] + nom_x[4];
   out_2012298189817593187[5] = delta_x[5] + nom_x[5];
   out_2012298189817593187[6] = delta_x[6] + nom_x[6];
   out_2012298189817593187[7] = delta_x[7] + nom_x[7];
   out_2012298189817593187[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6652027009148145501) {
   out_6652027009148145501[0] = -nom_x[0] + true_x[0];
   out_6652027009148145501[1] = -nom_x[1] + true_x[1];
   out_6652027009148145501[2] = -nom_x[2] + true_x[2];
   out_6652027009148145501[3] = -nom_x[3] + true_x[3];
   out_6652027009148145501[4] = -nom_x[4] + true_x[4];
   out_6652027009148145501[5] = -nom_x[5] + true_x[5];
   out_6652027009148145501[6] = -nom_x[6] + true_x[6];
   out_6652027009148145501[7] = -nom_x[7] + true_x[7];
   out_6652027009148145501[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6851681808913447295) {
   out_6851681808913447295[0] = 1.0;
   out_6851681808913447295[1] = 0;
   out_6851681808913447295[2] = 0;
   out_6851681808913447295[3] = 0;
   out_6851681808913447295[4] = 0;
   out_6851681808913447295[5] = 0;
   out_6851681808913447295[6] = 0;
   out_6851681808913447295[7] = 0;
   out_6851681808913447295[8] = 0;
   out_6851681808913447295[9] = 0;
   out_6851681808913447295[10] = 1.0;
   out_6851681808913447295[11] = 0;
   out_6851681808913447295[12] = 0;
   out_6851681808913447295[13] = 0;
   out_6851681808913447295[14] = 0;
   out_6851681808913447295[15] = 0;
   out_6851681808913447295[16] = 0;
   out_6851681808913447295[17] = 0;
   out_6851681808913447295[18] = 0;
   out_6851681808913447295[19] = 0;
   out_6851681808913447295[20] = 1.0;
   out_6851681808913447295[21] = 0;
   out_6851681808913447295[22] = 0;
   out_6851681808913447295[23] = 0;
   out_6851681808913447295[24] = 0;
   out_6851681808913447295[25] = 0;
   out_6851681808913447295[26] = 0;
   out_6851681808913447295[27] = 0;
   out_6851681808913447295[28] = 0;
   out_6851681808913447295[29] = 0;
   out_6851681808913447295[30] = 1.0;
   out_6851681808913447295[31] = 0;
   out_6851681808913447295[32] = 0;
   out_6851681808913447295[33] = 0;
   out_6851681808913447295[34] = 0;
   out_6851681808913447295[35] = 0;
   out_6851681808913447295[36] = 0;
   out_6851681808913447295[37] = 0;
   out_6851681808913447295[38] = 0;
   out_6851681808913447295[39] = 0;
   out_6851681808913447295[40] = 1.0;
   out_6851681808913447295[41] = 0;
   out_6851681808913447295[42] = 0;
   out_6851681808913447295[43] = 0;
   out_6851681808913447295[44] = 0;
   out_6851681808913447295[45] = 0;
   out_6851681808913447295[46] = 0;
   out_6851681808913447295[47] = 0;
   out_6851681808913447295[48] = 0;
   out_6851681808913447295[49] = 0;
   out_6851681808913447295[50] = 1.0;
   out_6851681808913447295[51] = 0;
   out_6851681808913447295[52] = 0;
   out_6851681808913447295[53] = 0;
   out_6851681808913447295[54] = 0;
   out_6851681808913447295[55] = 0;
   out_6851681808913447295[56] = 0;
   out_6851681808913447295[57] = 0;
   out_6851681808913447295[58] = 0;
   out_6851681808913447295[59] = 0;
   out_6851681808913447295[60] = 1.0;
   out_6851681808913447295[61] = 0;
   out_6851681808913447295[62] = 0;
   out_6851681808913447295[63] = 0;
   out_6851681808913447295[64] = 0;
   out_6851681808913447295[65] = 0;
   out_6851681808913447295[66] = 0;
   out_6851681808913447295[67] = 0;
   out_6851681808913447295[68] = 0;
   out_6851681808913447295[69] = 0;
   out_6851681808913447295[70] = 1.0;
   out_6851681808913447295[71] = 0;
   out_6851681808913447295[72] = 0;
   out_6851681808913447295[73] = 0;
   out_6851681808913447295[74] = 0;
   out_6851681808913447295[75] = 0;
   out_6851681808913447295[76] = 0;
   out_6851681808913447295[77] = 0;
   out_6851681808913447295[78] = 0;
   out_6851681808913447295[79] = 0;
   out_6851681808913447295[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2518437364328843337) {
   out_2518437364328843337[0] = state[0];
   out_2518437364328843337[1] = state[1];
   out_2518437364328843337[2] = state[2];
   out_2518437364328843337[3] = state[3];
   out_2518437364328843337[4] = state[4];
   out_2518437364328843337[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2518437364328843337[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2518437364328843337[7] = state[7];
   out_2518437364328843337[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6061066950155854042) {
   out_6061066950155854042[0] = 1;
   out_6061066950155854042[1] = 0;
   out_6061066950155854042[2] = 0;
   out_6061066950155854042[3] = 0;
   out_6061066950155854042[4] = 0;
   out_6061066950155854042[5] = 0;
   out_6061066950155854042[6] = 0;
   out_6061066950155854042[7] = 0;
   out_6061066950155854042[8] = 0;
   out_6061066950155854042[9] = 0;
   out_6061066950155854042[10] = 1;
   out_6061066950155854042[11] = 0;
   out_6061066950155854042[12] = 0;
   out_6061066950155854042[13] = 0;
   out_6061066950155854042[14] = 0;
   out_6061066950155854042[15] = 0;
   out_6061066950155854042[16] = 0;
   out_6061066950155854042[17] = 0;
   out_6061066950155854042[18] = 0;
   out_6061066950155854042[19] = 0;
   out_6061066950155854042[20] = 1;
   out_6061066950155854042[21] = 0;
   out_6061066950155854042[22] = 0;
   out_6061066950155854042[23] = 0;
   out_6061066950155854042[24] = 0;
   out_6061066950155854042[25] = 0;
   out_6061066950155854042[26] = 0;
   out_6061066950155854042[27] = 0;
   out_6061066950155854042[28] = 0;
   out_6061066950155854042[29] = 0;
   out_6061066950155854042[30] = 1;
   out_6061066950155854042[31] = 0;
   out_6061066950155854042[32] = 0;
   out_6061066950155854042[33] = 0;
   out_6061066950155854042[34] = 0;
   out_6061066950155854042[35] = 0;
   out_6061066950155854042[36] = 0;
   out_6061066950155854042[37] = 0;
   out_6061066950155854042[38] = 0;
   out_6061066950155854042[39] = 0;
   out_6061066950155854042[40] = 1;
   out_6061066950155854042[41] = 0;
   out_6061066950155854042[42] = 0;
   out_6061066950155854042[43] = 0;
   out_6061066950155854042[44] = 0;
   out_6061066950155854042[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6061066950155854042[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6061066950155854042[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6061066950155854042[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6061066950155854042[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6061066950155854042[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6061066950155854042[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6061066950155854042[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6061066950155854042[53] = -9.8000000000000007*dt;
   out_6061066950155854042[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6061066950155854042[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6061066950155854042[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6061066950155854042[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6061066950155854042[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6061066950155854042[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6061066950155854042[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6061066950155854042[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6061066950155854042[62] = 0;
   out_6061066950155854042[63] = 0;
   out_6061066950155854042[64] = 0;
   out_6061066950155854042[65] = 0;
   out_6061066950155854042[66] = 0;
   out_6061066950155854042[67] = 0;
   out_6061066950155854042[68] = 0;
   out_6061066950155854042[69] = 0;
   out_6061066950155854042[70] = 1;
   out_6061066950155854042[71] = 0;
   out_6061066950155854042[72] = 0;
   out_6061066950155854042[73] = 0;
   out_6061066950155854042[74] = 0;
   out_6061066950155854042[75] = 0;
   out_6061066950155854042[76] = 0;
   out_6061066950155854042[77] = 0;
   out_6061066950155854042[78] = 0;
   out_6061066950155854042[79] = 0;
   out_6061066950155854042[80] = 1;
}
void h_25(double *state, double *unused, double *out_2384782074905661785) {
   out_2384782074905661785[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3170119256419037021) {
   out_3170119256419037021[0] = 0;
   out_3170119256419037021[1] = 0;
   out_3170119256419037021[2] = 0;
   out_3170119256419037021[3] = 0;
   out_3170119256419037021[4] = 0;
   out_3170119256419037021[5] = 0;
   out_3170119256419037021[6] = 1;
   out_3170119256419037021[7] = 0;
   out_3170119256419037021[8] = 0;
}
void h_24(double *state, double *unused, double *out_4656752188981984794) {
   out_4656752188981984794[0] = state[4];
   out_4656752188981984794[1] = state[5];
}
void H_24(double *state, double *unused, double *out_387806342383286115) {
   out_387806342383286115[0] = 0;
   out_387806342383286115[1] = 0;
   out_387806342383286115[2] = 0;
   out_387806342383286115[3] = 0;
   out_387806342383286115[4] = 1;
   out_387806342383286115[5] = 0;
   out_387806342383286115[6] = 0;
   out_387806342383286115[7] = 0;
   out_387806342383286115[8] = 0;
   out_387806342383286115[9] = 0;
   out_387806342383286115[10] = 0;
   out_387806342383286115[11] = 0;
   out_387806342383286115[12] = 0;
   out_387806342383286115[13] = 0;
   out_387806342383286115[14] = 1;
   out_387806342383286115[15] = 0;
   out_387806342383286115[16] = 0;
   out_387806342383286115[17] = 0;
}
void h_30(double *state, double *unused, double *out_4203045851570633378) {
   out_4203045851570633378[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7697815586546645219) {
   out_7697815586546645219[0] = 0;
   out_7697815586546645219[1] = 0;
   out_7697815586546645219[2] = 0;
   out_7697815586546645219[3] = 0;
   out_7697815586546645219[4] = 1;
   out_7697815586546645219[5] = 0;
   out_7697815586546645219[6] = 0;
   out_7697815586546645219[7] = 0;
   out_7697815586546645219[8] = 0;
}
void h_26(double *state, double *unused, double *out_798627416090212267) {
   out_798627416090212267[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6911622575293093245) {
   out_6911622575293093245[0] = 0;
   out_6911622575293093245[1] = 0;
   out_6911622575293093245[2] = 0;
   out_6911622575293093245[3] = 0;
   out_6911622575293093245[4] = 0;
   out_6911622575293093245[5] = 0;
   out_6911622575293093245[6] = 0;
   out_6911622575293093245[7] = 1;
   out_6911622575293093245[8] = 0;
}
void h_27(double *state, double *unused, double *out_2616891192755183860) {
   out_2616891192755183860[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5474221515362702002) {
   out_5474221515362702002[0] = 0;
   out_5474221515362702002[1] = 0;
   out_5474221515362702002[2] = 0;
   out_5474221515362702002[3] = 1;
   out_5474221515362702002[4] = 0;
   out_5474221515362702002[5] = 0;
   out_5474221515362702002[6] = 0;
   out_5474221515362702002[7] = 0;
   out_5474221515362702002[8] = 0;
}
void h_29(double *state, double *unused, double *out_9204719119231479023) {
   out_9204719119231479023[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7187584242232253035) {
   out_7187584242232253035[0] = 0;
   out_7187584242232253035[1] = 1;
   out_7187584242232253035[2] = 0;
   out_7187584242232253035[3] = 0;
   out_7187584242232253035[4] = 0;
   out_7187584242232253035[5] = 0;
   out_7187584242232253035[6] = 0;
   out_7187584242232253035[7] = 0;
   out_7187584242232253035[8] = 0;
}
void h_28(double *state, double *unused, double *out_2724976043130846703) {
   out_2724976043130846703[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6176760814407768007) {
   out_6176760814407768007[0] = 1;
   out_6176760814407768007[1] = 0;
   out_6176760814407768007[2] = 0;
   out_6176760814407768007[3] = 0;
   out_6176760814407768007[4] = 0;
   out_6176760814407768007[5] = 0;
   out_6176760814407768007[6] = 0;
   out_6176760814407768007[7] = 0;
   out_6176760814407768007[8] = 0;
}
void h_31(double *state, double *unused, double *out_6022231603159684135) {
   out_6022231603159684135[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3139473294542076593) {
   out_3139473294542076593[0] = 0;
   out_3139473294542076593[1] = 0;
   out_3139473294542076593[2] = 0;
   out_3139473294542076593[3] = 0;
   out_3139473294542076593[4] = 0;
   out_3139473294542076593[5] = 0;
   out_3139473294542076593[6] = 0;
   out_3139473294542076593[7] = 0;
   out_3139473294542076593[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_2012298189817593187) {
  err_fun(nom_x, delta_x, out_2012298189817593187);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6652027009148145501) {
  inv_err_fun(nom_x, true_x, out_6652027009148145501);
}
void car_H_mod_fun(double *state, double *out_6851681808913447295) {
  H_mod_fun(state, out_6851681808913447295);
}
void car_f_fun(double *state, double dt, double *out_2518437364328843337) {
  f_fun(state,  dt, out_2518437364328843337);
}
void car_F_fun(double *state, double dt, double *out_6061066950155854042) {
  F_fun(state,  dt, out_6061066950155854042);
}
void car_h_25(double *state, double *unused, double *out_2384782074905661785) {
  h_25(state, unused, out_2384782074905661785);
}
void car_H_25(double *state, double *unused, double *out_3170119256419037021) {
  H_25(state, unused, out_3170119256419037021);
}
void car_h_24(double *state, double *unused, double *out_4656752188981984794) {
  h_24(state, unused, out_4656752188981984794);
}
void car_H_24(double *state, double *unused, double *out_387806342383286115) {
  H_24(state, unused, out_387806342383286115);
}
void car_h_30(double *state, double *unused, double *out_4203045851570633378) {
  h_30(state, unused, out_4203045851570633378);
}
void car_H_30(double *state, double *unused, double *out_7697815586546645219) {
  H_30(state, unused, out_7697815586546645219);
}
void car_h_26(double *state, double *unused, double *out_798627416090212267) {
  h_26(state, unused, out_798627416090212267);
}
void car_H_26(double *state, double *unused, double *out_6911622575293093245) {
  H_26(state, unused, out_6911622575293093245);
}
void car_h_27(double *state, double *unused, double *out_2616891192755183860) {
  h_27(state, unused, out_2616891192755183860);
}
void car_H_27(double *state, double *unused, double *out_5474221515362702002) {
  H_27(state, unused, out_5474221515362702002);
}
void car_h_29(double *state, double *unused, double *out_9204719119231479023) {
  h_29(state, unused, out_9204719119231479023);
}
void car_H_29(double *state, double *unused, double *out_7187584242232253035) {
  H_29(state, unused, out_7187584242232253035);
}
void car_h_28(double *state, double *unused, double *out_2724976043130846703) {
  h_28(state, unused, out_2724976043130846703);
}
void car_H_28(double *state, double *unused, double *out_6176760814407768007) {
  H_28(state, unused, out_6176760814407768007);
}
void car_h_31(double *state, double *unused, double *out_6022231603159684135) {
  h_31(state, unused, out_6022231603159684135);
}
void car_H_31(double *state, double *unused, double *out_3139473294542076593) {
  H_31(state, unused, out_3139473294542076593);
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
