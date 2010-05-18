#include <math.h>
#include <stdio.h>
 
#include "macros_and_parameters.h"
 
 
/*
 
   Fitting Formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies.
   Daniel J. Eisenstein & Wayne Hu, Institute for Advanced Study
 
   There are two primary routines here, one to set the cosmology, the
   other to construct the transfer function for a single wavenumber k.
   You should call the former once (per cosmology) and the latter as
   many times as you want.
 
   TFmdm_set_cosm() -- User passes all the cosmological parameters as
	arguments; the routine sets up all of the scalar quantites needed
	computation of the fitting formula.  The input parameters are:
	1) omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density.
	2) omega_baryon -- Density of baryons, in units of critical.
	3) omega_hdm    -- Density of massive neutrinos, in units of critical
	4) degen_hdm    -- (Int) Number of degenerate massive neutrino species
	5) omega_lambda -- Cosmological constant
	6) hubble       -- Hubble constant, in units of 100 km/s/Mpc
	7) redshift     -- The redshift at which to evaluate
 
   TFmdm_onek_mpc() -- User passes a single wavenumber, in units of Mpc^-1.
	Routine returns the transfer function from the Eisenstein & Hu
	fitting formula, based on the cosmology currently held in the
	internal variables.  The routine returns T_cb (the CDM+Baryon
	density-weighted transfer function), although T_cbn (the CDM+
	Baryon+Neutrino density-weighted transfer function) is stored
	in the global variable tf_cbnu.
 
   We also supply TFmdm_onek_hmpc(), which is identical to the previous
	routine, but takes the wavenumber in units of h Mpc^-1.
 
   We hold the internal scalar quantities in global variables, so that
   the user may access them in an external program, via "extern" declarations.
 
   Please note that all internal length scales are in Mpc, not h^-1 Mpc!
 
*/
 
/* -------------------------- Prototypes ----------------------------- */
 
 
int TFmdm_set_cosm(FLOAT omega_matter, FLOAT omega_baryon, FLOAT omega_hdm,
	int degen_hdm, FLOAT omega_lambda, FLOAT hubble, FLOAT redshift);
FLOAT TFmdm_onek_mpc(FLOAT kk);
FLOAT TFmdm_onek_hmpc(FLOAT kk);
FLOAT TF98_onek_mpc(FLOAT kk);
 
// Convenience from Numerical Recipes in C, 2nd edition
static FLOAT sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
 
/* ------------------------- Global Variables ------------------------ */
 
// The following are set in TFmdm_set_cosm()
 
FLOAT   alpha_gamma,	/* sqrt(alpha_nu) */
	alpha_nu,	/* The small-scale suppression */
	beta_c,		/* The correction to the log in the small-scale */
	num_degen_hdm,	/* Number of degenerate massive neutrino species */
	f_baryon,	/* Baryon fraction */
	f_bnu,		/* Baryon + Massive Neutrino fraction */
	f_cb,		/* Baryon + CDM fraction */
	f_cdm,		/* CDM fraction */
	f_hdm,		/* Massive Neutrino fraction */
	growth_k0,	/* D_1(z) -- the growth function as k->0 */
	growth_to_z0,	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
	hhubble,	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
	k_equality,	/* The comoving wave number of the horizon at equality*/
        R_equality,	/* R_eq parameter */
        z_equality,	/* Redshift of matter-radiation equality */
	obhh,		/* Omega_baryon * hubble^2 */
	ochh,		/* Omega_cdm * hubble^2 */
	omega_curv,	/* = 1 - omega_matter - omega_lambda */
	omega_lambda_z, /* Omega_lambda at the given redshift */
	omega_matter_z,	/* Omega_matter at the given redshift */
	omhh,		/* Omega_matter * hubble^2 */
	onhh,		/* Omega_hdm * hubble^2 */
	p_c,		/* The correction to the exponent before drag epoch */
	p_cb,		/* The correction to the exponent after drag epoch */
	sound_horizon_fit,  /* The sound horizon at the drag epoch */
	sound_horizon,  /* The sound horizon at the drag epoch */
        k_Silk,         /* The Silk damping scale */
	theta_cmb,	/* The temperature of the CMB, in units of 2.7 K */
	y_drag,		/* Ratio of z_equality to z_drag */
	z_drag,		/* Redshift of the drag epoch */
        R_drag;	        /* R_drag parameter */
 
// The following are set in TFmdm_onek_mpc()
 
FLOAT	gamma_eff,	/* Effective \Gamma */
	growth_cb,	/* Growth factor for CDM+Baryon perturbations */
	growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
	max_fs_correction,  /* Correction near maximal free streaming */
	qq,		/* Wavenumber rescaled by \Gamma */
	qq_eff,		/* Wavenumber rescaled by effective Gamma */
	qq_nu,		/* Wavenumber compared to maximal free streaming */
	tf_master,	/* Master TF */
	tf_sup,		/* Suppressed TF */
	y_freestream; 	/* The epoch of free-streaming for a given scale */
 
// Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as
 
FLOAT   tf_cb,		/* The transfer function for density-weighted
			CDM + Baryon perturbations. */
	tf_cbnu;	/* The transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */
 
/* By default, these functions return tf_cb */
 
/* ------------------------- TFmdm_set_cosm() ------------------------ */
 
int TFmdm_set_cosm(FLOAT omega_matter, FLOAT omega_baryon, FLOAT omega_hdm,
	int degen_hdm, FLOAT omega_lambda, FLOAT hubble, FLOAT redshift)
 
/* This routine takes cosmological parameters and a redshift and sets up
all the internal scalar quantities needed to compute the transfer function. */
 
/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density. */
/* 	  omega_baryon -- Density of baryons, in units of critical. */
/* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
/* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
/*        omega_lambda -- Cosmological constant */
/* 	  hubble       -- Hubble constant, in units of 100 km/s/Mpc */
/*        redshift     -- The redshift at which to evaluate */
 
/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
	sets many global variables for use in TFmdm_onek_mpc() */
{
    FLOAT z_drag_b1, z_drag_b2, omega_denom;
    int qwarn;
    qwarn = 0;
 
    theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */
 
    /* Look for strange input */
    if (omega_baryon<0.0) {
	fprintf(stderr,
	  "TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
	qwarn = 1;
    }
    if (omega_hdm<0.0) {
	fprintf(stderr,
	  "TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n");
	qwarn = 1;
    }
    if (hubble<=0.0) {
	fprintf(stderr,"TFmdm_set_cosm(): Negative Hubble constant illegal.\n");
	return(1);  /* Can't recover */
    } else if (hubble>2.0) {
	fprintf(stderr,"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
	qwarn = 1;
    }
    if (redshift<=-1.0) {
	fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
	return(1);
    } else if (redshift>99.0) {
      //	fprintf(stderr,
		// "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
	//	qwarn = 1;
    }
    if (degen_hdm<1) degen_hdm=1;
    num_degen_hdm = (FLOAT) degen_hdm;	
	/* Have to save this for TFmdm_onek_mpc() */
    /* This routine would crash if baryons or neutrinos were zero,
	so don't allow that */
    if (omega_baryon<=0) omega_baryon=1e-5;
    if (omega_hdm<=0) omega_hdm=1e-5;
 
    omega_curv = 1.0-omega_matter-omega_lambda;
    omhh = omega_matter*SQR(hubble);
    obhh = omega_baryon*SQR(hubble);
    ochh = (omega_matter - omega_baryon)*SQR(hubble);
    onhh = omega_hdm*SQR(hubble);
    f_baryon = omega_baryon/omega_matter;
    f_hdm = omega_hdm/omega_matter;
    f_cdm = 1.0-f_baryon-f_hdm;
    f_cb = f_cdm+f_baryon;
    f_bnu = f_baryon+f_hdm;
 
    /* Compute the equality scale. */
    z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
    k_equality = 0.0746*omhh/SQR(theta_cmb);
    R_equality = 31.5 * obhh / SQR(SQR(theta_cmb)) * (1.0e3 / z_equality);

    /* Compute the drag epoch and sound horizon */
    z_drag_b1 = 0.313*POW(omhh,-0.419)*(1+0.607*POW(omhh,0.674));
    z_drag_b2 = 0.238*POW(omhh,0.223);
    z_drag = 1291*POW(omhh,0.251)/(1.0+0.659*POW(omhh,0.828))*
		(1.0+z_drag_b1*POW(obhh,z_drag_b2));
    R_drag = 31.5 * obhh / SQR(SQR(theta_cmb)) * (1.0e3 / z_drag);
    y_drag = z_equality/(1.0+z_drag);

    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*POW(obhh,0.75));
 
    sound_horizon = 2.0/(3.0*k_equality) * sqrt(6.0/R_equality) * log( (sqrt(1.0+R_drag) + sqrt(R_drag + R_equality)) / (1.0 + sqrt(R_equality)) );
      

    /* Silk damping scale */
    k_Silk = 1.6 * POW(obhh,0.52) * POW(omhh,0.73) * ( 1.0 + POW(10.4*omhh,-0.95) );


    /* Set up for the free-streaming & infall growth function */
    p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));
 
    omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+
			omega_matter*(1.0+redshift));
    omega_lambda_z = omega_lambda/omega_denom;
    omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
    growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/
	    (POW(omega_matter_z,4.0/7.0)-omega_lambda_z+
	    (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
    growth_to_z0 = z_equality*2.5*omega_matter/(POW(omega_matter,4.0/7.0)
	    -omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
    growth_to_z0 = growth_k0/growth_to_z0;	
 
    /* Compute small-scale suppression */
    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
	POW(1+y_drag,p_cb-p_c)*
	(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
	(1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*POW(num_degen_hdm,0.2))*
	(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = sqrt(alpha_nu);
    beta_c = 1/(1-0.949*f_bnu);
    /* Done setting scalar variables */
    hhubble = hubble;	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
    return qwarn;
}
 
/* ---------------------------- TFmdm_onek_mpc() ---------------------- */
 
FLOAT TFmdm_onek_mpc(FLOAT kk)
 
/* Given a wavenumber in Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
 
/* Input: kk -- Wavenumber in Mpc^-1 */
 
/* Output: The following are set as global variables:
	growth_cb -- the transfer function for density-weighted
			CDM + Baryon perturbations.
 	growth_cbnu -- the transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */
{
    FLOAT tf_sup_L, tf_sup_C;
    FLOAT temp1, temp2;
 
    qq = kk/omhh*SQR(theta_cmb);
 
    /* Compute the scale-dependent growth functions */
    y_freestream = 17.2*f_hdm*(1+0.488*POW(f_hdm,-7.0/6.0))*
		SQR(num_degen_hdm*qq/f_hdm);
    temp1 = POW(growth_k0, 1.0-p_cb);
    temp2 = POW(growth_k0/(1+y_freestream),0.7);
    growth_cb = POW(1.0+temp2, p_cb/0.7)*temp1;
    growth_cbnu = POW(POW(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;
 
    /* Compute the master function */
    gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
		(1+SQR(SQR(kk*sound_horizon_fit*0.43))));
    qq_eff = qq*omhh/gamma_eff;
 
    tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
    tf_sup_C = 14.4+325/(1+60.5*POW(qq_eff,1.11));
    tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));
 
    qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
    max_fs_correction = 1+1.2*POW(f_hdm,0.64)*POW(num_degen_hdm,0.3+0.6*f_hdm)/
		(POW(qq_nu,-1.6)+POW(qq_nu,0.8));
    tf_master = tf_sup*max_fs_correction;
 
    /* Now compute the CDM+HDM+baryon transfer functions */
    tf_cb = tf_master*growth_cb/growth_k0;
    tf_cbnu = tf_master*growth_cbnu/growth_k0;
    return tf_cb;
}
 
/* ---------------------------- TFmdm_onek_hmpc() ---------------------- */
 
FLOAT TFmdm_onek_hmpc(FLOAT kk)
 
/* Given a wavenumber in h Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
 
/* Input: kk -- Wavenumber in h Mpc^-1 */
 
/* Output: The following are set as global variables:
	growth_cb -- the transfer function for density-weighted
			CDM + Baryon perturbations.
 	growth_cbnu -- the transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */
{
    return TFmdm_onek_mpc(kk*hhubble);
}



/* ---------------------------- TF98_onek_mpc() ---------------------- */
/*
  This is the version from the 1998 Eisenstein & Hu paper
  (1998ApJ...496..605E), which includes the baryonic wiggles.  This
  transfer function can be selected by setting PowerSpectrumType=15.

  Added by Michael Kuhlen, February 2010.
*/ 
FLOAT TF98_onek_mpc(FLOAT kk)
 
/* Given a wavenumber in Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
 
/* Input: kk -- Wavenumber in Mpc^-1 */
 
/* The function returns transfer_cb */
{
  FLOAT a1, a2, alpha_c, b1, b2, beta_c;
  FLOAT q, G, alpha_b, ks, f;
  FLOAT C1, Tc1, C2, Tc2, Tc;
  FLOAT beta_node, stilde, beta_b, Tb1, Tb;
  
  q = kk / (13.41 * k_equality);

  a1 = POW(46.9 * omhh,0.670) * ( 1.0 + POW(32.1 * omhh,-0.532) );
  a2 = POW(12.0 * omhh,0.424) * ( 1.0 + (45.0 * omhh,-0.582) );
  alpha_c = POW(a1,-obhh/omhh) * POW(a2,-POW(obhh/omhh,3));
  
  b1 = 0.944 / ( 1.0 + POW(458.0 * omhh,-0.708) );
  b2 = POW(0.395 * omhh,-0.0266);
  beta_c = 1.0 / ( 1.0 + b1 * ( POW(ochh/omhh,b2) - 1.0 ) );

  G = y_drag * ( -6.0 * sqrt(1.0 + y_drag) + (2.0 + 3.0*y_drag) * log( (sqrt(1.0+y_drag) + 1.0) / (sqrt(1.0+y_drag) - 1.0) ) );

  alpha_b = 2.07 * k_equality * sound_horizon * POW(1.0 + R_drag,-0.75) * G;


  ks = kk * sound_horizon;
  f = 1.0 / (1.0 + POW(ks/5.4,4));
  
  C1 = 14.2 + 386.0 / (1.0 + 69.9*POW(q,1.08));
  Tc1 = log(exp(1.0) + 1.8 * beta_c * q) / (log(exp(1.0) + 1.8 * beta_c * q) + C1*q*q);
  
  C2 = 14.2/alpha_c + 386.0 / (1.0 + 69.9*POW(q,1.08));
  Tc2 = log(exp(1.0) + 1.8 * beta_c *q) / (log(exp(1.0) + 1.8 * beta_c * q) + C2*q*q);

  Tc = f * Tc1 + (1.0 - f) * Tc2;


  beta_node = 8.41 * POW(omhh,0.435);
  stilde = sound_horizon / POW(1.0 + POW(beta_node / ks,3),1.0/3.0);

  beta_b = 0.5 + obhh/omhh + (3.0 - 2.0*obhh/omhh)*sqrt(POW(17.2*omhh,2) + 1.0);

  Tb1 = log(exp(1.0) + 1.8 * q) / (log(exp(1.0) + 1.8 * q) + C1*q*q);

  Tb = ( Tb1 / (1.0 + POW(ks/5.2,2)) + alpha_b / (1.0 + POW(beta_b/ks,3)) * exp(-POW(kk/k_Silk,1.4)) ) * sin(kk*stilde) / (kk*stilde);

  
  tf_cb = obhh/omhh * Tb + ochh/omhh * Tc;
  
  return tf_cb;
}
