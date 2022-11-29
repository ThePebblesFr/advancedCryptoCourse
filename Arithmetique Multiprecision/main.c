#include "elliptic_curves.h"

int main()
{
	printf("--------------------------- BEGINNING PROCESS ---------------------------\n\n");
	printf("----------------------------- ELLIPTIC CURVE ----------------------------\n");
	Curve* C = malloc(sizeof(Curve));

	mpz_t a, b, p;
	mpz_inits(a, b, p, NULL);

	mpz_set_ui(a, 0);
	mpz_set_ui(b, 3); // Y²=X³+3 avec p = 7
	mpz_set_ui(p, 7);

	gmp_printf("---> Initializing Curve C : y² = x³ + %Zdx + %Zd, avec (x,y) appartenant à F%Zd\n", a, b, p);
	initializeCurve(C, a, b, p);


	printf("-------------------- CARTESIAN ADDITION AND DOUBLING --------------------\n");

	P_carthesian* P = malloc(sizeof(P_carthesian));
	P_carthesian* Q = malloc(sizeof(P_carthesian));
	P_carthesian* R = malloc(sizeof(P_carthesian));
	mpz_t x1, x2, y1, y2, z1, z2, k;
	mpz_inits(x1, x2, y1, y2, z1, z2, k, NULL);

	mpz_set_ui(x1, 6);
	mpz_set_ui(y1, 4); // P = (1, 2)
	mpz_set_ui(x2, 1); // Q = (4, 5)
	mpz_set_ui(y2, 5);

	gmp_printf("---> Initializing P = (%Zd, %Zd)\n", x1, y1);
	initializeP_Cartesian(P, x1, y1, C);
	gmp_printf("---> Initializing Q = (%Zd, %Zd)\n", x2, y2);
	initializeP_Cartesian(Q, x2, y2, C);

	printf("---> Computing P + Q\n");
	addition(R, P, Q, C);
	gmp_printf("[CARTESIAN-]   P + Q = (%Zd, %Zd)\n", R->x, R->y);

	printf("---> Computing 2P\n");
	doubling(R, P, C);
	gmp_printf("[CARTESIAN-]      2P = (%Zd, %Zd)\n", R->x, R->y);

	printf("--------------------- JACOBIAN ADDITION AND DOUBLING --------------------\n");
	P_jacobian* Pj = malloc(sizeof(P_jacobian));
	P_jacobian* Qj = malloc(sizeof(P_jacobian));
	P_jacobian* Rj = malloc(sizeof(P_jacobian));
	mpz_set_ui(x1, 5);
	mpz_set_ui(y1, 3); // Pj = (5, 3, 4)
	mpz_set_ui(z1, 4);
	mpz_set_ui(x2, 1);
	mpz_set_ui(y2, 2); // Qj = (1, 2, 1)
	mpz_set_ui(z2, 1);

	gmp_printf("---> Initializing Pj = (%Zd, %Zd, %Zd)\n", x1, y1, z1);
	initializeP_Jacobian(Pj, x1, y1, z1, C);
	gmp_printf("---> Initializing Qj = (%Zd, %Zd, %Zd)\n", x2, y2, z2);
	initializeP_Jacobian(Qj, x2, y2, z2, C);

	printf("---> Computing Pj + Qj\n");
	addition_jacob(Rj, Pj, Qj, C);
	gmp_printf("[JACOBIAN--] Pj + Qj = (%Zd, %Zd, %Zd)\n", Rj->X, Rj->Y, Rj->Z);

	printf("---> Computing 2Pj\n");
	doubling_jacob(Rj, Pj, C);
	gmp_printf("[JACOBIAN--]     2Pj = (%Zd, %Zd, %Zd)\n", Rj->X, Rj->Y, Rj->Z);

	printf("---> Reinitializing P\n");
	mpz_set_ui(x1, 1); // P = (1, 5)
	mpz_set_ui(y1, 5);
	mpz_set_ui(k, 2);  // k = 2
	initializeP_Cartesian(P, x1, y1, C);

	printf("-------------------------- SCALAR MULTIPLICATION ------------------------\n");
	gmp_printf("---> Computing kP = %Zd(%Zd, %Zd) using Montgomery scale\n", k, P->x, P->y);
	montgomeryScale(R, P, k, C);
	gmp_printf("[CARTESIAN-]      %ZdP = (%Zd, %Zd)\n", k, R->x, R->y);

	printf("\n----------------------------- ENDING PROCESS ----------------------------\n");
	mpz_clears(a, b, p, x1, x2, y1, y2, z1, z2, k, NULL);
	return 0;
}