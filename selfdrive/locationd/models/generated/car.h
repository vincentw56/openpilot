#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_8738564341327856279);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_858109794346859647);
void car_H_mod_fun(double *state, double *out_3333624422478311790);
void car_f_fun(double *state, double dt, double *out_4907349842632933819);
void car_F_fun(double *state, double dt, double *out_8981105644751692590);
void car_h_25(double *state, double *unused, double *out_5988176754801101980);
void car_H_25(double *state, double *unused, double *out_8759459417119999305);
void car_h_24(double *state, double *unused, double *out_821681974249550378);
void car_H_24(double *state, double *unused, double *out_6533751633141130743);
void car_h_30(double *state, double *unused, double *out_6263370817085607869);
void car_H_30(double *state, double *unused, double *out_1842769075628382550);
void car_h_26(double *state, double *unused, double *out_8867776137894305520);
void car_H_26(double *state, double *unused, double *out_8102605353009687401);
void car_h_27(double *state, double *unused, double *out_410864762612747166);
void car_H_27(double *state, double *unused, double *out_4017532387428807461);
void car_h_29(double *state, double *unused, double *out_3450818713815754486);
void car_H_29(double *state, double *unused, double *out_1332537731313990366);
void car_h_28(double *state, double *unused, double *out_1354374615514409601);
void car_H_28(double *state, double *unused, double *out_6414936748383520940);
void car_h_31(double *state, double *unused, double *out_3961906637644709644);
void car_H_31(double *state, double *unused, double *out_8728813455243038877);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}