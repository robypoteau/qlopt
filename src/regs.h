#ifndef REGS_H
#define REGS_H

#include <dbg.h>
#include <misc.h>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <numerical_integration.h>
#include <mp_spline.h>

mp_mat mpFindA(const mp_vec& t, const mp_mat& U, int m);
mpreal mpPsuedoDet(const mp_mat& A);
mp_vec mpFindP(const mp_vec& t, const mp_mat& U, const mp_vec& dx, int m);
mpreal mpInnerProd(const mp_vec& u1, const mp_vec& u2, const mp_vec& time);
vec regularization(soln_env *env);

mp_mat mpQLinearRungeKutta4(string fname, const mp_vec& time, const mp_vec& u, const mp_vec& yNot, const mp_mat& xNminus);
mpreal mpSimpson(const mp_vec& t, const mp_vec& x);
mp_mat mpQlinear(mp_sys fhandle, const mpreal& t, const mp_vec& x, const mp_vec& u, thesis::mp_spline* Xn, int n);

mpreal mp_norm(const mp_mat& M);
mp_mat mp_inverse(const mp_mat& M);
mpreal mp_cond(const mp_mat& A);
mp_mat mpReshape(const mp_mat& U, int n, int m);

mp_mat mpCofactor(const mp_mat& A);
mp_mat mpRowColRemoval(const mp_mat& A, const int row, const int col);
#endif