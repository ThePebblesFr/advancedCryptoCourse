#include "elliptic_curves.h"

void initializeCurve(Curve* c, mpz_t a, mpz_t b, mpz_t p)
{
	mpz_set(c->a, a);
	mpz_set(c->b, b);
	mpz_set(c->p, p);
}

void inverseModulaire(mpz_t x_inv, mpz_t x, Curve* c)
{
	// a^(p-2)%p
	mpz_t exposant, res_puissance;
	mpz_inits(exposant, res_puissance, NULL);
	mpz_sub_ui(exposant, c->p, 2);
	mpz_powm(x_inv, x, exposant, c->p);

	mpz_clears(exposant, res_puissance, NULL);
}

// void fromCartesianToJacobian(Polynome* P, Curve* c)
// {
// 	mpz_t X_Jacobian, Y_Jacobian;
// 	mpz_inits(X_Jacobian, Y_Jacobian, NULL);
// 	mpz_mul_ui(X_Jacobian, P->cartesian_coord[0], pow(CHOSEN_Z, 2));
// 	mpz_mul_ui(Y_Jacobian, P->cartesian_coord[1], pow(CHOSEN_Z, 3));
	
// 	mpz_mod(P->jacobian_coord[0], X_Jacobian, c->p);
// 	mpz_mod(P->jacobian_coord[1], Y_Jacobian, c->p);
// 	mpz_set_ui(P->jacobian_coord[2], CHOSEN_Z);

// 	mpz_clears(X_Jacobian, Y_Jacobian, NULL);
// }

// void fromJacobianToCartesian(Polynome* P, Curve* c)
// {
// 	mpz_t X_Cartesian, Y_Cartesian, InvZ, InvZ_2, InvZ_3;
// 	mpz_inits(X_Cartesian, Y_Cartesian, InvZ, InvZ_2, InvZ_3, NULL);
// 	inverseModulaire(InvZ, P->jacobian_coord[2], c);

// 	mpz_powm_ui(InvZ_2, InvZ, 2, c->p);
// 	mpz_powm_ui(InvZ_3, InvZ, 3, c->p);

// 	mpz_mul(X_Cartesian, P->jacobian_coord[0], InvZ_2);
// 	mpz_mul(Y_Cartesian, P->jacobian_coord[1], InvZ_3);

// 	mpz_mod(P->cartesian_coord[0], X_Cartesian, c->p);
// 	mpz_mod(P->cartesian_coord[1], X_Cartesian, c->p);

// 	mpz_clears(X_Cartesian, Y_Cartesian, InvZ, InvZ_2, InvZ_3, NULL);
// }

void initializeP_Cartesian(P_carthesian* P, mpz_t x, mpz_t y, Curve* c)
{
	mpz_inits(P->x, P->y, NULL);
	mpz_mod(P->x, x, c->p);
	mpz_mod(P->y, y, c->p);
}

void addition(P_carthesian* R, P_carthesian* P, P_carthesian* Q, Curve* c)
{
	mpz_t lambda, X2_X1, Y2_Y1, Inv_X2_X1;
	mpz_inits(lambda, X2_X1, Y2_Y1, Inv_X2_X1, NULL);

	// Computing lambda
	mpz_sub(X2_X1, Q->x, P->x);
	mpz_sub(Y2_Y1, Q->y, P->y);
	inverseModulaire(Inv_X2_X1, X2_X1, c);
	mpz_mul(lambda, Y2_Y1, Inv_X2_X1);
	mpz_mod(lambda, lambda, c->p);

	// Computing X3
	mpz_powm_ui(R->x, lambda, 2, c->p);
	mpz_sub(R->x, R->x, P->x);
	mpz_sub(R->x, R->x, Q->x);
	mpz_mod(R->x, R->x, c->p);

	// Computing Y3
	mpz_sub(R->y, P->x, R->x);
	mpz_mul(R->y, lambda, R->y);
	mpz_sub(R->y, R->y, P->y);
	mpz_mod(R->y, R->y, c->p);

	mpz_clears(lambda, X2_X1, Y2_Y1, Inv_X2_X1, NULL);
}

void doubling(P_carthesian* R, P_carthesian* P, Curve* c)
{
	mpz_t lambda, inter1, Inv_2Y1;
	mpz_inits(lambda, inter1, Inv_2Y1, NULL);

	// Computing lambda
	mpz_mul_ui(inter1, P->y, 2);
	inverseModulaire(Inv_2Y1, inter1, c);
	mpz_powm_ui(inter1, P->x, 2, c->p);
	mpz_mul_ui(inter1, inter1, 3);
	mpz_add(inter1, inter1, c->a);
	mpz_mul(lambda, inter1, Inv_2Y1);
	mpz_mod(lambda, lambda, c->p);

	// Computing x3
	mpz_powm_ui(R->x, lambda, 2, c->p);
	mpz_mul_ui(inter1, P->x, 2);
	mpz_sub(R->x, R->x, inter1);
	mpz_mod(R->x, R->x, c->p);

	// Computing y3
	mpz_sub(R->y, P->x, R->x);
	mpz_mul(R->y, lambda, R->y);
	mpz_sub(R->y, R->y, P->y);
	mpz_mod(R->y, R->y, c->p);

	mpz_clears(lambda, inter1, Inv_2Y1, NULL);
}

void initializeP_Jacobian(P_jacobian* P, mpz_t X, mpz_t Y, mpz_t Z, Curve* c)
{
	mpz_inits(P->X, P->Y, P->Z, NULL);
	mpz_mod(P->X, X, c->p);
	mpz_mod(P->Y, Y, c->p);
	mpz_mod(P->Z, Z, c->p);
}

void addition_jacob(P_jacobian* R, P_jacobian* P, P_jacobian* Q, Curve* c)
{
	mpz_t Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, r, inter1, V;
	mpz_inits(Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, r, inter1, V, NULL);

	// Compute Z1Z1
	mpz_powm_ui(Z1Z1, P->Z, 2, c->p);
	// Compute Z2Z2
	mpz_powm_ui(Z2Z2, Q->Z, 2, c->p);
	// Compute U1
	mpz_mul(U1, P->X, Z2Z2);
	mpz_mod(U1, U1, c->p);
	// Compute U2
	mpz_mul(U2, Q->X, Z1Z1);
	mpz_mod(U2, U2, c->p);
	// Compute S1
	mpz_mul(S1, P->Y, Q->Z);
	mpz_mod(S1, S1, c->p);
	mpz_mul(S1, S1, Z2Z2);
	mpz_mod(S1, S1, c->p);
	// Compute S2
	mpz_mul(S2, Q->Y, P->Z);
	mpz_mod(S2, S2, c->p);
	mpz_mul(S2, S2, Z1Z1);
	mpz_mod(S2, S2, c->p);
	// Compute H
	mpz_sub(H, U2, U1);
	mpz_mod(H, H, c->p);
	// Compute I
	mpz_mul_ui(I, H, 2);
	mpz_mod(I, I, c->p);
	mpz_powm_ui(I, I, 2, c->p);
	// Compute J
	mpz_mul(J, H, I);
	mpz_mod(J, J, c->p);
	// Compute r
	mpz_sub(r, S2, S1);
	mpz_mod(r, r, c->p);
	mpz_mul_ui(r, r, 2);
	mpz_mod(r, r, c->p);
	// Compute V
	mpz_mul(V, U1, I);
	mpz_mod(V, V, c->p);
	// Compute X3
	mpz_powm_ui(R->X, r, 2, c->p);
	mpz_mod(R->X, R->X, c->p);
	mpz_sub(R->X, R->X, J);
	mpz_mod(R->X, R->X, c->p);
	mpz_mul_ui(inter1, V, 2);
	mpz_mod(inter1, inter1, c->p);
	mpz_sub(R->X, R->X, inter1);
	mpz_mod(R->X, R->X, c->p);
	// Compute Y3
	mpz_sub(R->Y, V, R->X);
	mpz_mod(R->Y, R->Y, c->p);
	mpz_mul(R->Y, r, R->Y);
	mpz_mod(R->Y, R->Y, c->p);
	mpz_mul(inter1, S1, J);
	mpz_mod(inter1, inter1, c->p);
	mpz_mul_ui(inter1, inter1, 2);
	mpz_mod(inter1, inter1, c->p);
	mpz_sub(R->Y, R->Y, inter1);
	mpz_mod(R->Y, R->Y, c->p);
	// Compute Z3
	mpz_add(R->Z, P->Z, Q->Z);
	mpz_mod(R->Z, R->Z, c->p);
	mpz_powm_ui(R->Z, R->Z, 2, c->p);
	mpz_sub(R->Z, R->Z, Z1Z1);
	mpz_mod(R->Z, R->Z, c->p);
	mpz_sub(R->Z, R->Z, Z2Z2);
	mpz_mod(R->Z, R->Z, c->p);
	mpz_mul(R->Z, R->Z, H);
	mpz_mod(R->Z, R->Z, c->p);

	mpz_clears(Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, r, inter1, V, NULL);
}

void doubling_jacob(P_jacobian* R, P_jacobian* P, Curve* c)
{
	mpz_t A, B, C, D, E, F, inter1;
	mpz_inits(A, B, C, D, E, F, inter1, NULL);

	// Compute A
	mpz_powm_ui(A, P->X, 2, c->p);
	// Compute B
	mpz_powm_ui(B, P->Y, 2, c->p);
	// Compute C
	mpz_powm_ui(C, B, 2, c->p);
	// Compute D
	mpz_add(D, P->X, B);
	mpz_mod(D, D, c->p);
	mpz_powm_ui(D, D, 2, c->p);
	mpz_sub(D, D, A);
	mpz_mod(D, D, c->p);
	mpz_sub(D, D, C);
	mpz_mod(D, D, c->p);
	mpz_mul_ui(D, D, 2);
	mpz_mod(D, D, c->p);
	// Compute E
	mpz_mul_ui(E, A, 3);
	mpz_mod(E, E, c->p);
	// Compute F
	mpz_powm_ui(F, E, 2, c->p);
	// Compute X3
	mpz_mul_ui(R->X, D, 2);
	mpz_mod(R->X, R->X, c->p);
	mpz_sub(R->X, F, R->X);
	mpz_mod(R->X, R->X, c->p);
	// Compute Y3
	mpz_sub(R->Y, D, R->X);
	mpz_mod(R->Y, R->Y, c->p);
	mpz_mul(R->Y, E, R->Y);
	mpz_mod(R->Y, R->Y, c->p);
	mpz_mul_ui(inter1, C, 8);
	mpz_mod(inter1, inter1, c->p);
	mpz_sub(R->Y, R->Y, inter1);
	mpz_mod(R->Y, R->Y, c->p);
	// Compute Z3
	mpz_mul(R->Z, P->Y, P->Z);
	mpz_mod(R->Z, R->Z, c->p);
	mpz_mul_ui(R->Z, R->Z, 2);
	mpz_mod(R->Z, R->Z, c->p);
}

void montgomeryScale(P_carthesian* R, P_carthesian* P, mpz_t k, Curve* c)
{
	int n = mpz_sizeinbase(k, 2);
	P_carthesian returned_P[2], R_inter1, R_inter2;
	mpz_inits(returned_P[0].x, returned_P[0].y, returned_P[1].x, returned_P[1].y, R_inter1.x, R_inter1.y, R_inter2.x, R_inter2.y, NULL);

	mpz_set(returned_P[0].x, P->x);
	mpz_set(returned_P[0].y, P->y);
	doubling(&returned_P[1], P, c);
	int b = 0;

	for (int i = n - 2; i >= 0; i--)
	{
		b = mpz_tstbit(k, i);
		addition(&R_inter1, returned_P + 1 - b, returned_P + b, c);
		doubling(&R_inter2, returned_P + b, c);

		mpz_set(returned_P[1 - b].x, R_inter1.x);
		mpz_set(returned_P[1 - b].y, R_inter1.y);
		mpz_set(returned_P[b].x, R_inter2.x);
		mpz_set(returned_P[b].y, R_inter2.y);
	}
	mpz_set(R->x, returned_P[0].x);
	mpz_set(R->y, returned_P[0].y);

	mpz_clears(returned_P[0].x, returned_P[0].y, returned_P[1].x, returned_P[1].y, R_inter1.x, R_inter1.y, R_inter2.x, R_inter2.y, NULL);
}
