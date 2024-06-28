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
void err_fun(double *nom_x, double *delta_x, double *out_6748780409524464664) {
   out_6748780409524464664[0] = delta_x[0] + nom_x[0];
   out_6748780409524464664[1] = delta_x[1] + nom_x[1];
   out_6748780409524464664[2] = delta_x[2] + nom_x[2];
   out_6748780409524464664[3] = delta_x[3] + nom_x[3];
   out_6748780409524464664[4] = delta_x[4] + nom_x[4];
   out_6748780409524464664[5] = delta_x[5] + nom_x[5];
   out_6748780409524464664[6] = delta_x[6] + nom_x[6];
   out_6748780409524464664[7] = delta_x[7] + nom_x[7];
   out_6748780409524464664[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7365381153033288851) {
   out_7365381153033288851[0] = -nom_x[0] + true_x[0];
   out_7365381153033288851[1] = -nom_x[1] + true_x[1];
   out_7365381153033288851[2] = -nom_x[2] + true_x[2];
   out_7365381153033288851[3] = -nom_x[3] + true_x[3];
   out_7365381153033288851[4] = -nom_x[4] + true_x[4];
   out_7365381153033288851[5] = -nom_x[5] + true_x[5];
   out_7365381153033288851[6] = -nom_x[6] + true_x[6];
   out_7365381153033288851[7] = -nom_x[7] + true_x[7];
   out_7365381153033288851[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4852154337395239742) {
   out_4852154337395239742[0] = 1.0;
   out_4852154337395239742[1] = 0;
   out_4852154337395239742[2] = 0;
   out_4852154337395239742[3] = 0;
   out_4852154337395239742[4] = 0;
   out_4852154337395239742[5] = 0;
   out_4852154337395239742[6] = 0;
   out_4852154337395239742[7] = 0;
   out_4852154337395239742[8] = 0;
   out_4852154337395239742[9] = 0;
   out_4852154337395239742[10] = 1.0;
   out_4852154337395239742[11] = 0;
   out_4852154337395239742[12] = 0;
   out_4852154337395239742[13] = 0;
   out_4852154337395239742[14] = 0;
   out_4852154337395239742[15] = 0;
   out_4852154337395239742[16] = 0;
   out_4852154337395239742[17] = 0;
   out_4852154337395239742[18] = 0;
   out_4852154337395239742[19] = 0;
   out_4852154337395239742[20] = 1.0;
   out_4852154337395239742[21] = 0;
   out_4852154337395239742[22] = 0;
   out_4852154337395239742[23] = 0;
   out_4852154337395239742[24] = 0;
   out_4852154337395239742[25] = 0;
   out_4852154337395239742[26] = 0;
   out_4852154337395239742[27] = 0;
   out_4852154337395239742[28] = 0;
   out_4852154337395239742[29] = 0;
   out_4852154337395239742[30] = 1.0;
   out_4852154337395239742[31] = 0;
   out_4852154337395239742[32] = 0;
   out_4852154337395239742[33] = 0;
   out_4852154337395239742[34] = 0;
   out_4852154337395239742[35] = 0;
   out_4852154337395239742[36] = 0;
   out_4852154337395239742[37] = 0;
   out_4852154337395239742[38] = 0;
   out_4852154337395239742[39] = 0;
   out_4852154337395239742[40] = 1.0;
   out_4852154337395239742[41] = 0;
   out_4852154337395239742[42] = 0;
   out_4852154337395239742[43] = 0;
   out_4852154337395239742[44] = 0;
   out_4852154337395239742[45] = 0;
   out_4852154337395239742[46] = 0;
   out_4852154337395239742[47] = 0;
   out_4852154337395239742[48] = 0;
   out_4852154337395239742[49] = 0;
   out_4852154337395239742[50] = 1.0;
   out_4852154337395239742[51] = 0;
   out_4852154337395239742[52] = 0;
   out_4852154337395239742[53] = 0;
   out_4852154337395239742[54] = 0;
   out_4852154337395239742[55] = 0;
   out_4852154337395239742[56] = 0;
   out_4852154337395239742[57] = 0;
   out_4852154337395239742[58] = 0;
   out_4852154337395239742[59] = 0;
   out_4852154337395239742[60] = 1.0;
   out_4852154337395239742[61] = 0;
   out_4852154337395239742[62] = 0;
   out_4852154337395239742[63] = 0;
   out_4852154337395239742[64] = 0;
   out_4852154337395239742[65] = 0;
   out_4852154337395239742[66] = 0;
   out_4852154337395239742[67] = 0;
   out_4852154337395239742[68] = 0;
   out_4852154337395239742[69] = 0;
   out_4852154337395239742[70] = 1.0;
   out_4852154337395239742[71] = 0;
   out_4852154337395239742[72] = 0;
   out_4852154337395239742[73] = 0;
   out_4852154337395239742[74] = 0;
   out_4852154337395239742[75] = 0;
   out_4852154337395239742[76] = 0;
   out_4852154337395239742[77] = 0;
   out_4852154337395239742[78] = 0;
   out_4852154337395239742[79] = 0;
   out_4852154337395239742[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8529546029908871000) {
   out_8529546029908871000[0] = state[0];
   out_8529546029908871000[1] = state[1];
   out_8529546029908871000[2] = state[2];
   out_8529546029908871000[3] = state[3];
   out_8529546029908871000[4] = state[4];
   out_8529546029908871000[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8529546029908871000[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8529546029908871000[7] = state[7];
   out_8529546029908871000[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1356409578556224397) {
   out_1356409578556224397[0] = 1;
   out_1356409578556224397[1] = 0;
   out_1356409578556224397[2] = 0;
   out_1356409578556224397[3] = 0;
   out_1356409578556224397[4] = 0;
   out_1356409578556224397[5] = 0;
   out_1356409578556224397[6] = 0;
   out_1356409578556224397[7] = 0;
   out_1356409578556224397[8] = 0;
   out_1356409578556224397[9] = 0;
   out_1356409578556224397[10] = 1;
   out_1356409578556224397[11] = 0;
   out_1356409578556224397[12] = 0;
   out_1356409578556224397[13] = 0;
   out_1356409578556224397[14] = 0;
   out_1356409578556224397[15] = 0;
   out_1356409578556224397[16] = 0;
   out_1356409578556224397[17] = 0;
   out_1356409578556224397[18] = 0;
   out_1356409578556224397[19] = 0;
   out_1356409578556224397[20] = 1;
   out_1356409578556224397[21] = 0;
   out_1356409578556224397[22] = 0;
   out_1356409578556224397[23] = 0;
   out_1356409578556224397[24] = 0;
   out_1356409578556224397[25] = 0;
   out_1356409578556224397[26] = 0;
   out_1356409578556224397[27] = 0;
   out_1356409578556224397[28] = 0;
   out_1356409578556224397[29] = 0;
   out_1356409578556224397[30] = 1;
   out_1356409578556224397[31] = 0;
   out_1356409578556224397[32] = 0;
   out_1356409578556224397[33] = 0;
   out_1356409578556224397[34] = 0;
   out_1356409578556224397[35] = 0;
   out_1356409578556224397[36] = 0;
   out_1356409578556224397[37] = 0;
   out_1356409578556224397[38] = 0;
   out_1356409578556224397[39] = 0;
   out_1356409578556224397[40] = 1;
   out_1356409578556224397[41] = 0;
   out_1356409578556224397[42] = 0;
   out_1356409578556224397[43] = 0;
   out_1356409578556224397[44] = 0;
   out_1356409578556224397[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1356409578556224397[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1356409578556224397[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1356409578556224397[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1356409578556224397[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1356409578556224397[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1356409578556224397[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1356409578556224397[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1356409578556224397[53] = -9.8000000000000007*dt;
   out_1356409578556224397[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1356409578556224397[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1356409578556224397[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1356409578556224397[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1356409578556224397[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1356409578556224397[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1356409578556224397[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1356409578556224397[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1356409578556224397[62] = 0;
   out_1356409578556224397[63] = 0;
   out_1356409578556224397[64] = 0;
   out_1356409578556224397[65] = 0;
   out_1356409578556224397[66] = 0;
   out_1356409578556224397[67] = 0;
   out_1356409578556224397[68] = 0;
   out_1356409578556224397[69] = 0;
   out_1356409578556224397[70] = 1;
   out_1356409578556224397[71] = 0;
   out_1356409578556224397[72] = 0;
   out_1356409578556224397[73] = 0;
   out_1356409578556224397[74] = 0;
   out_1356409578556224397[75] = 0;
   out_1356409578556224397[76] = 0;
   out_1356409578556224397[77] = 0;
   out_1356409578556224397[78] = 0;
   out_1356409578556224397[79] = 0;
   out_1356409578556224397[80] = 1;
}
void h_25(double *state, double *unused, double *out_7811054891928739484) {
   out_7811054891928739484[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8060467551967462791) {
   out_8060467551967462791[0] = 0;
   out_8060467551967462791[1] = 0;
   out_8060467551967462791[2] = 0;
   out_8060467551967462791[3] = 0;
   out_8060467551967462791[4] = 0;
   out_8060467551967462791[5] = 0;
   out_8060467551967462791[6] = 1;
   out_8060467551967462791[7] = 0;
   out_8060467551967462791[8] = 0;
}
void h_24(double *state, double *unused, double *out_7092083788774129740) {
   out_7092083788774129740[0] = state[4];
   out_7092083788774129740[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3240146047311474528) {
   out_3240146047311474528[0] = 0;
   out_3240146047311474528[1] = 0;
   out_3240146047311474528[2] = 0;
   out_3240146047311474528[3] = 0;
   out_3240146047311474528[4] = 1;
   out_3240146047311474528[5] = 0;
   out_3240146047311474528[6] = 0;
   out_3240146047311474528[7] = 0;
   out_3240146047311474528[8] = 0;
   out_3240146047311474528[9] = 0;
   out_3240146047311474528[10] = 0;
   out_3240146047311474528[11] = 0;
   out_3240146047311474528[12] = 0;
   out_3240146047311474528[13] = 0;
   out_3240146047311474528[14] = 1;
   out_3240146047311474528[15] = 0;
   out_3240146047311474528[16] = 0;
   out_3240146047311474528[17] = 0;
}
void h_30(double *state, double *unused, double *out_47789635108067043) {
   out_47789635108067043[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7931128604824222721) {
   out_7931128604824222721[0] = 0;
   out_7931128604824222721[1] = 0;
   out_7931128604824222721[2] = 0;
   out_7931128604824222721[3] = 0;
   out_7931128604824222721[4] = 1;
   out_7931128604824222721[5] = 0;
   out_7931128604824222721[6] = 0;
   out_7931128604824222721[7] = 0;
   out_7931128604824222721[8] = 0;
}
void h_26(double *state, double *unused, double *out_1656372417640759219) {
   out_1656372417640759219[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4318964233093406567) {
   out_4318964233093406567[0] = 0;
   out_4318964233093406567[1] = 0;
   out_4318964233093406567[2] = 0;
   out_4318964233093406567[3] = 0;
   out_4318964233093406567[4] = 0;
   out_4318964233093406567[5] = 0;
   out_4318964233093406567[6] = 0;
   out_4318964233093406567[7] = 1;
   out_4318964233093406567[8] = 0;
}
void h_27(double *state, double *unused, double *out_2831809747985136497) {
   out_2831809747985136497[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5756365293023797810) {
   out_5756365293023797810[0] = 0;
   out_5756365293023797810[1] = 0;
   out_5756365293023797810[2] = 0;
   out_5756365293023797810[3] = 1;
   out_5756365293023797810[4] = 0;
   out_5756365293023797810[5] = 0;
   out_5756365293023797810[6] = 0;
   out_5756365293023797810[7] = 0;
   out_5756365293023797810[8] = 0;
}
void h_29(double *state, double *unused, double *out_5049462902768912688) {
   out_5049462902768912688[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8441359949138614905) {
   out_8441359949138614905[0] = 0;
   out_8441359949138614905[1] = 1;
   out_8441359949138614905[2] = 0;
   out_8441359949138614905[3] = 0;
   out_8441359949138614905[4] = 0;
   out_8441359949138614905[5] = 0;
   out_8441359949138614905[6] = 0;
   out_8441359949138614905[7] = 0;
   out_8441359949138614905[8] = 0;
}
void h_28(double *state, double *unused, double *out_2169863519675709148) {
   out_2169863519675709148[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1039396450915283797) {
   out_1039396450915283797[0] = 1;
   out_1039396450915283797[1] = 0;
   out_1039396450915283797[2] = 0;
   out_1039396450915283797[3] = 0;
   out_1039396450915283797[4] = 0;
   out_1039396450915283797[5] = 0;
   out_1039396450915283797[6] = 0;
   out_1039396450915283797[7] = 0;
   out_1039396450915283797[8] = 0;
}
void h_31(double *state, double *unused, double *out_70217758825309701) {
   out_70217758825309701[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8091113513844423219) {
   out_8091113513844423219[0] = 0;
   out_8091113513844423219[1] = 0;
   out_8091113513844423219[2] = 0;
   out_8091113513844423219[3] = 0;
   out_8091113513844423219[4] = 0;
   out_8091113513844423219[5] = 0;
   out_8091113513844423219[6] = 0;
   out_8091113513844423219[7] = 0;
   out_8091113513844423219[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_6748780409524464664) {
  err_fun(nom_x, delta_x, out_6748780409524464664);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7365381153033288851) {
  inv_err_fun(nom_x, true_x, out_7365381153033288851);
}
void car_H_mod_fun(double *state, double *out_4852154337395239742) {
  H_mod_fun(state, out_4852154337395239742);
}
void car_f_fun(double *state, double dt, double *out_8529546029908871000) {
  f_fun(state,  dt, out_8529546029908871000);
}
void car_F_fun(double *state, double dt, double *out_1356409578556224397) {
  F_fun(state,  dt, out_1356409578556224397);
}
void car_h_25(double *state, double *unused, double *out_7811054891928739484) {
  h_25(state, unused, out_7811054891928739484);
}
void car_H_25(double *state, double *unused, double *out_8060467551967462791) {
  H_25(state, unused, out_8060467551967462791);
}
void car_h_24(double *state, double *unused, double *out_7092083788774129740) {
  h_24(state, unused, out_7092083788774129740);
}
void car_H_24(double *state, double *unused, double *out_3240146047311474528) {
  H_24(state, unused, out_3240146047311474528);
}
void car_h_30(double *state, double *unused, double *out_47789635108067043) {
  h_30(state, unused, out_47789635108067043);
}
void car_H_30(double *state, double *unused, double *out_7931128604824222721) {
  H_30(state, unused, out_7931128604824222721);
}
void car_h_26(double *state, double *unused, double *out_1656372417640759219) {
  h_26(state, unused, out_1656372417640759219);
}
void car_H_26(double *state, double *unused, double *out_4318964233093406567) {
  H_26(state, unused, out_4318964233093406567);
}
void car_h_27(double *state, double *unused, double *out_2831809747985136497) {
  h_27(state, unused, out_2831809747985136497);
}
void car_H_27(double *state, double *unused, double *out_5756365293023797810) {
  H_27(state, unused, out_5756365293023797810);
}
void car_h_29(double *state, double *unused, double *out_5049462902768912688) {
  h_29(state, unused, out_5049462902768912688);
}
void car_H_29(double *state, double *unused, double *out_8441359949138614905) {
  H_29(state, unused, out_8441359949138614905);
}
void car_h_28(double *state, double *unused, double *out_2169863519675709148) {
  h_28(state, unused, out_2169863519675709148);
}
void car_H_28(double *state, double *unused, double *out_1039396450915283797) {
  H_28(state, unused, out_1039396450915283797);
}
void car_h_31(double *state, double *unused, double *out_70217758825309701) {
  h_31(state, unused, out_70217758825309701);
}
void car_H_31(double *state, double *unused, double *out_8091113513844423219) {
  H_31(state, unused, out_8091113513844423219);
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
