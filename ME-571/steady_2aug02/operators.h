#ifndef _OPERATORS_
#define _OPERATORS_

void Eval_Au (double *u, double *A_on_u);
void Eval_RHS (double *u, double *rhs);
void To_Build (double *u, double *u_prime);

void B_op (double *u);
void Green (double *u, double *u_prime);
void Make_Operators (double *u_base, double omega_base, double dr);
void Operators_ini (void);
void Operators_fin (void);

#endif /* ifndef _OPERATORS_ */
