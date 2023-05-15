/******************************************
/
/ Chem Equilibrium Data Tables
/
/*****************************************/

struct EquilibriumTableType {

int dim_size;

// Arrays of size dim_size
double* density;
double* temperature;

// Arrays of size dim_size**2
double* HI;
double* HII;
double* DI;
double* DII;
double* HM;
double* H2I;
double* H2II;
double* HDI;
double* HeI;
double* HeII;
double* HeIII;
double* de;

};
