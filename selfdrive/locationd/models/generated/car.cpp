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
void err_fun(double *nom_x, double *delta_x, double *out_8738564341327856279) {
   out_8738564341327856279[0] = delta_x[0] + nom_x[0];
   out_8738564341327856279[1] = delta_x[1] + nom_x[1];
   out_8738564341327856279[2] = delta_x[2] + nom_x[2];
   out_8738564341327856279[3] = delta_x[3] + nom_x[3];
   out_8738564341327856279[4] = delta_x[4] + nom_x[4];
   out_8738564341327856279[5] = delta_x[5] + nom_x[5];
   out_8738564341327856279[6] = delta_x[6] + nom_x[6];
   out_8738564341327856279[7] = delta_x[7] + nom_x[7];
   out_8738564341327856279[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_858109794346859647) {
   out_858109794346859647[0] = -nom_x[0] + true_x[0];
   out_858109794346859647[1] = -nom_x[1] + true_x[1];
   out_858109794346859647[2] = -nom_x[2] + true_x[2];
   out_858109794346859647[3] = -nom_x[3] + true_x[3];
   out_858109794346859647[4] = -nom_x[4] + true_x[4];
   out_858109794346859647[5] = -nom_x[5] + true_x[5];
   out_858109794346859647[6] = -nom_x[6] + true_x[6];
   out_858109794346859647[7] = -nom_x[7] + true_x[7];
   out_858109794346859647[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3333624422478311790) {
   out_3333624422478311790[0] = 1.0;
   out_3333624422478311790[1] = 0;
   out_3333624422478311790[2] = 0;
   out_3333624422478311790[3] = 0;
   out_3333624422478311790[4] = 0;
   out_3333624422478311790[5] = 0;
   out_3333624422478311790[6] = 0;
   out_3333624422478311790[7] = 0;
   out_3333624422478311790[8] = 0;
   out_3333624422478311790[9] = 0;
   out_3333624422478311790[10] = 1.0;
   out_3333624422478311790[11] = 0;
   out_3333624422478311790[12] = 0;
   out_3333624422478311790[13] = 0;
   out_3333624422478311790[14] = 0;
   out_3333624422478311790[15] = 0;
   out_3333624422478311790[16] = 0;
   out_3333624422478311790[17] = 0;
   out_3333624422478311790[18] = 0;
   out_3333624422478311790[19] = 0;
   out_3333624422478311790[20] = 1.0;
   out_3333624422478311790[21] = 0;
   out_3333624422478311790[22] = 0;
   out_3333624422478311790[23] = 0;
   out_3333624422478311790[24] = 0;
   out_3333624422478311790[25] = 0;
   out_3333624422478311790[26] = 0;
   out_3333624422478311790[27] = 0;
   out_3333624422478311790[28] = 0;
   out_3333624422478311790[29] = 0;
   out_3333624422478311790[30] = 1.0;
   out_3333624422478311790[31] = 0;
   out_3333624422478311790[32] = 0;
   out_3333624422478311790[33] = 0;
   out_3333624422478311790[34] = 0;
   out_3333624422478311790[35] = 0;
   out_3333624422478311790[36] = 0;
   out_3333624422478311790[37] = 0;
   out_3333624422478311790[38] = 0;
   out_3333624422478311790[39] = 0;
   out_3333624422478311790[40] = 1.0;
   out_3333624422478311790[41] = 0;
   out_3333624422478311790[42] = 0;
   out_3333624422478311790[43] = 0;
   out_3333624422478311790[44] = 0;
   out_3333624422478311790[45] = 0;
   out_3333624422478311790[46] = 0;
   out_3333624422478311790[47] = 0;
   out_3333624422478311790[48] = 0;
   out_3333624422478311790[49] = 0;
   out_3333624422478311790[50] = 1.0;
   out_3333624422478311790[51] = 0;
   out_3333624422478311790[52] = 0;
   out_3333624422478311790[53] = 0;
   out_3333624422478311790[54] = 0;
   out_3333624422478311790[55] = 0;
   out_3333624422478311790[56] = 0;
   out_3333624422478311790[57] = 0;
   out_3333624422478311790[58] = 0;
   out_3333624422478311790[59] = 0;
   out_3333624422478311790[60] = 1.0;
   out_3333624422478311790[61] = 0;
   out_3333624422478311790[62] = 0;
   out_3333624422478311790[63] = 0;
   out_3333624422478311790[64] = 0;
   out_3333624422478311790[65] = 0;
   out_3333624422478311790[66] = 0;
   out_3333624422478311790[67] = 0;
   out_3333624422478311790[68] = 0;
   out_3333624422478311790[69] = 0;
   out_3333624422478311790[70] = 1.0;
   out_3333624422478311790[71] = 0;
   out_3333624422478311790[72] = 0;
   out_3333624422478311790[73] = 0;
   out_3333624422478311790[74] = 0;
   out_3333624422478311790[75] = 0;
   out_3333624422478311790[76] = 0;
   out_3333624422478311790[77] = 0;
   out_3333624422478311790[78] = 0;
   out_3333624422478311790[79] = 0;
   out_3333624422478311790[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4907349842632933819) {
   out_4907349842632933819[0] = state[0];
   out_4907349842632933819[1] = state[1];
   out_4907349842632933819[2] = state[2];
   out_4907349842632933819[3] = state[3];
   out_4907349842632933819[4] = state[4];
   out_4907349842632933819[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4907349842632933819[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4907349842632933819[7] = state[7];
   out_4907349842632933819[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8981105644751692590) {
   out_8981105644751692590[0] = 1;
   out_8981105644751692590[1] = 0;
   out_8981105644751692590[2] = 0;
   out_8981105644751692590[3] = 0;
   out_8981105644751692590[4] = 0;
   out_8981105644751692590[5] = 0;
   out_8981105644751692590[6] = 0;
   out_8981105644751692590[7] = 0;
   out_8981105644751692590[8] = 0;
   out_8981105644751692590[9] = 0;
   out_8981105644751692590[10] = 1;
   out_8981105644751692590[11] = 0;
   out_8981105644751692590[12] = 0;
   out_8981105644751692590[13] = 0;
   out_8981105644751692590[14] = 0;
   out_8981105644751692590[15] = 0;
   out_8981105644751692590[16] = 0;
   out_8981105644751692590[17] = 0;
   out_8981105644751692590[18] = 0;
   out_8981105644751692590[19] = 0;
   out_8981105644751692590[20] = 1;
   out_8981105644751692590[21] = 0;
   out_8981105644751692590[22] = 0;
   out_8981105644751692590[23] = 0;
   out_8981105644751692590[24] = 0;
   out_8981105644751692590[25] = 0;
   out_8981105644751692590[26] = 0;
   out_8981105644751692590[27] = 0;
   out_8981105644751692590[28] = 0;
   out_8981105644751692590[29] = 0;
   out_8981105644751692590[30] = 1;
   out_8981105644751692590[31] = 0;
   out_8981105644751692590[32] = 0;
   out_8981105644751692590[33] = 0;
   out_8981105644751692590[34] = 0;
   out_8981105644751692590[35] = 0;
   out_8981105644751692590[36] = 0;
   out_8981105644751692590[37] = 0;
   out_8981105644751692590[38] = 0;
   out_8981105644751692590[39] = 0;
   out_8981105644751692590[40] = 1;
   out_8981105644751692590[41] = 0;
   out_8981105644751692590[42] = 0;
   out_8981105644751692590[43] = 0;
   out_8981105644751692590[44] = 0;
   out_8981105644751692590[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8981105644751692590[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8981105644751692590[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8981105644751692590[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8981105644751692590[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8981105644751692590[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8981105644751692590[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8981105644751692590[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8981105644751692590[53] = -9.8000000000000007*dt;
   out_8981105644751692590[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8981105644751692590[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8981105644751692590[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8981105644751692590[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8981105644751692590[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8981105644751692590[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8981105644751692590[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8981105644751692590[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8981105644751692590[62] = 0;
   out_8981105644751692590[63] = 0;
   out_8981105644751692590[64] = 0;
   out_8981105644751692590[65] = 0;
   out_8981105644751692590[66] = 0;
   out_8981105644751692590[67] = 0;
   out_8981105644751692590[68] = 0;
   out_8981105644751692590[69] = 0;
   out_8981105644751692590[70] = 1;
   out_8981105644751692590[71] = 0;
   out_8981105644751692590[72] = 0;
   out_8981105644751692590[73] = 0;
   out_8981105644751692590[74] = 0;
   out_8981105644751692590[75] = 0;
   out_8981105644751692590[76] = 0;
   out_8981105644751692590[77] = 0;
   out_8981105644751692590[78] = 0;
   out_8981105644751692590[79] = 0;
   out_8981105644751692590[80] = 1;
}
void h_25(double *state, double *unused, double *out_5988176754801101980) {
   out_5988176754801101980[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8759459417119999305) {
   out_8759459417119999305[0] = 0;
   out_8759459417119999305[1] = 0;
   out_8759459417119999305[2] = 0;
   out_8759459417119999305[3] = 0;
   out_8759459417119999305[4] = 0;
   out_8759459417119999305[5] = 0;
   out_8759459417119999305[6] = 1;
   out_8759459417119999305[7] = 0;
   out_8759459417119999305[8] = 0;
}
void h_24(double *state, double *unused, double *out_821681974249550378) {
   out_821681974249550378[0] = state[4];
   out_821681974249550378[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6533751633141130743) {
   out_6533751633141130743[0] = 0;
   out_6533751633141130743[1] = 0;
   out_6533751633141130743[2] = 0;
   out_6533751633141130743[3] = 0;
   out_6533751633141130743[4] = 1;
   out_6533751633141130743[5] = 0;
   out_6533751633141130743[6] = 0;
   out_6533751633141130743[7] = 0;
   out_6533751633141130743[8] = 0;
   out_6533751633141130743[9] = 0;
   out_6533751633141130743[10] = 0;
   out_6533751633141130743[11] = 0;
   out_6533751633141130743[12] = 0;
   out_6533751633141130743[13] = 0;
   out_6533751633141130743[14] = 1;
   out_6533751633141130743[15] = 0;
   out_6533751633141130743[16] = 0;
   out_6533751633141130743[17] = 0;
}
void h_30(double *state, double *unused, double *out_6263370817085607869) {
   out_6263370817085607869[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1842769075628382550) {
   out_1842769075628382550[0] = 0;
   out_1842769075628382550[1] = 0;
   out_1842769075628382550[2] = 0;
   out_1842769075628382550[3] = 0;
   out_1842769075628382550[4] = 1;
   out_1842769075628382550[5] = 0;
   out_1842769075628382550[6] = 0;
   out_1842769075628382550[7] = 0;
   out_1842769075628382550[8] = 0;
}
void h_26(double *state, double *unused, double *out_8867776137894305520) {
   out_8867776137894305520[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8102605353009687401) {
   out_8102605353009687401[0] = 0;
   out_8102605353009687401[1] = 0;
   out_8102605353009687401[2] = 0;
   out_8102605353009687401[3] = 0;
   out_8102605353009687401[4] = 0;
   out_8102605353009687401[5] = 0;
   out_8102605353009687401[6] = 0;
   out_8102605353009687401[7] = 1;
   out_8102605353009687401[8] = 0;
}
void h_27(double *state, double *unused, double *out_410864762612747166) {
   out_410864762612747166[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4017532387428807461) {
   out_4017532387428807461[0] = 0;
   out_4017532387428807461[1] = 0;
   out_4017532387428807461[2] = 0;
   out_4017532387428807461[3] = 1;
   out_4017532387428807461[4] = 0;
   out_4017532387428807461[5] = 0;
   out_4017532387428807461[6] = 0;
   out_4017532387428807461[7] = 0;
   out_4017532387428807461[8] = 0;
}
void h_29(double *state, double *unused, double *out_3450818713815754486) {
   out_3450818713815754486[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1332537731313990366) {
   out_1332537731313990366[0] = 0;
   out_1332537731313990366[1] = 1;
   out_1332537731313990366[2] = 0;
   out_1332537731313990366[3] = 0;
   out_1332537731313990366[4] = 0;
   out_1332537731313990366[5] = 0;
   out_1332537731313990366[6] = 0;
   out_1332537731313990366[7] = 0;
   out_1332537731313990366[8] = 0;
}
void h_28(double *state, double *unused, double *out_1354374615514409601) {
   out_1354374615514409601[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6414936748383520940) {
   out_6414936748383520940[0] = 1;
   out_6414936748383520940[1] = 0;
   out_6414936748383520940[2] = 0;
   out_6414936748383520940[3] = 0;
   out_6414936748383520940[4] = 0;
   out_6414936748383520940[5] = 0;
   out_6414936748383520940[6] = 0;
   out_6414936748383520940[7] = 0;
   out_6414936748383520940[8] = 0;
}
void h_31(double *state, double *unused, double *out_3961906637644709644) {
   out_3961906637644709644[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8728813455243038877) {
   out_8728813455243038877[0] = 0;
   out_8728813455243038877[1] = 0;
   out_8728813455243038877[2] = 0;
   out_8728813455243038877[3] = 0;
   out_8728813455243038877[4] = 0;
   out_8728813455243038877[5] = 0;
   out_8728813455243038877[6] = 0;
   out_8728813455243038877[7] = 0;
   out_8728813455243038877[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_8738564341327856279) {
  err_fun(nom_x, delta_x, out_8738564341327856279);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_858109794346859647) {
  inv_err_fun(nom_x, true_x, out_858109794346859647);
}
void car_H_mod_fun(double *state, double *out_3333624422478311790) {
  H_mod_fun(state, out_3333624422478311790);
}
void car_f_fun(double *state, double dt, double *out_4907349842632933819) {
  f_fun(state,  dt, out_4907349842632933819);
}
void car_F_fun(double *state, double dt, double *out_8981105644751692590) {
  F_fun(state,  dt, out_8981105644751692590);
}
void car_h_25(double *state, double *unused, double *out_5988176754801101980) {
  h_25(state, unused, out_5988176754801101980);
}
void car_H_25(double *state, double *unused, double *out_8759459417119999305) {
  H_25(state, unused, out_8759459417119999305);
}
void car_h_24(double *state, double *unused, double *out_821681974249550378) {
  h_24(state, unused, out_821681974249550378);
}
void car_H_24(double *state, double *unused, double *out_6533751633141130743) {
  H_24(state, unused, out_6533751633141130743);
}
void car_h_30(double *state, double *unused, double *out_6263370817085607869) {
  h_30(state, unused, out_6263370817085607869);
}
void car_H_30(double *state, double *unused, double *out_1842769075628382550) {
  H_30(state, unused, out_1842769075628382550);
}
void car_h_26(double *state, double *unused, double *out_8867776137894305520) {
  h_26(state, unused, out_8867776137894305520);
}
void car_H_26(double *state, double *unused, double *out_8102605353009687401) {
  H_26(state, unused, out_8102605353009687401);
}
void car_h_27(double *state, double *unused, double *out_410864762612747166) {
  h_27(state, unused, out_410864762612747166);
}
void car_H_27(double *state, double *unused, double *out_4017532387428807461) {
  H_27(state, unused, out_4017532387428807461);
}
void car_h_29(double *state, double *unused, double *out_3450818713815754486) {
  h_29(state, unused, out_3450818713815754486);
}
void car_H_29(double *state, double *unused, double *out_1332537731313990366) {
  H_29(state, unused, out_1332537731313990366);
}
void car_h_28(double *state, double *unused, double *out_1354374615514409601) {
  h_28(state, unused, out_1354374615514409601);
}
void car_H_28(double *state, double *unused, double *out_6414936748383520940) {
  H_28(state, unused, out_6414936748383520940);
}
void car_h_31(double *state, double *unused, double *out_3961906637644709644) {
  h_31(state, unused, out_3961906637644709644);
}
void car_H_31(double *state, double *unused, double *out_8728813455243038877) {
  H_31(state, unused, out_8728813455243038877);
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
