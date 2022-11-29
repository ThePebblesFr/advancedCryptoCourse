#ifndef ELLIPTIC_CURVES_H_
#define ELLIPTIC_CURVES_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

#define CHOSEN_Z 1 // value for z when passing from cartesian to jacobian coordinates

typedef struct {
	mpz_t a;
	mpz_t b;
	mpz_t p;
} Curve;

typedef struct {
	mpz_t x;
	mpz_t y;
} P_carthesian;

typedef struct {
	mpz_t X;
	mpz_t Y;
	mpz_t Z;
} P_jacobian;

void initializeCurve(Curve* c, mpz_t a, mpz_t b, mpz_t p);
void inverseModulaire(mpz_t x_inv, mpz_t x, Curve* c);
// void fromCartesianToJacobian(Polynome* P, Curve* c);
// void fromJacobianToCartesian(Polynome* P, Curve* c);

void initializeP_Cartesian(P_carthesian* P, mpz_t x, mpz_t y, Curve* c);
void addition(P_carthesian* R, P_carthesian* P, P_carthesian* Q, Curve* c);
void doubling(P_carthesian* R, P_carthesian* P, Curve* c);

void initializeP_Jacobian(P_jacobian* P, mpz_t X, mpz_t Y, mpz_t Z, Curve* c);
void addition_jacob(P_jacobian* R, P_jacobian* P, P_jacobian* Q, Curve* c);
void doubling_jacob(P_jacobian* R, P_jacobian* P, Curve* c);

void montgomeryScale(P_carthesian* R, P_carthesian* P,  mpz_t k, Curve* c);

#endif