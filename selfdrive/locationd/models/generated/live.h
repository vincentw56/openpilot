#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8634487822663781268);
void live_err_fun(double *nom_x, double *delta_x, double *out_783016460155138911);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_106834797478900715);
void live_H_mod_fun(double *state, double *out_6579876420749052012);
void live_f_fun(double *state, double dt, double *out_3659162376472982505);
void live_F_fun(double *state, double dt, double *out_4600865459334655732);
void live_h_4(double *state, double *unused, double *out_2796482808698657408);
void live_H_4(double *state, double *unused, double *out_5953790806247218602);
void live_h_9(double *state, double *unused, double *out_7894829364926935980);
void live_H_9(double *state, double *unused, double *out_6194980452876809247);
void live_h_10(double *state, double *unused, double *out_8392480068734059084);
void live_H_10(double *state, double *unused, double *out_8023184458012303454);
void live_h_12(double *state, double *unused, double *out_3668127217955138174);
void live_H_12(double *state, double *unused, double *out_6574889831294812269);
void live_h_35(double *state, double *unused, double *out_7475238483394447684);
void live_H_35(double *state, double *unused, double *out_9126291210089725638);
void live_h_32(double *state, double *unused, double *out_412463461389678539);
void live_H_32(double *state, double *unused, double *out_2405588359925709622);
void live_h_13(double *state, double *unused, double *out_2678666661056351279);
void live_H_13(double *state, double *unused, double *out_3784475594372073238);
void live_h_14(double *state, double *unused, double *out_7894829364926935980);
void live_H_14(double *state, double *unused, double *out_6194980452876809247);
void live_h_33(double *state, double *unused, double *out_1066723490886474754);
void live_H_33(double *state, double *unused, double *out_5975734205450868034);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}