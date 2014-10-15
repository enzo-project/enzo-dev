
/***********************************************************************
/
/  STOCHASTIC FORCING CLASS METHOD: Init
/
/  written by: Wolfram Schmidt
/  date:       May, 2005
/  modified1: Oct, 2014: updated to support Enzo 2.4 // P. Grete
			 using huge_number instead of FLT_MAX
/
/  PURPOSE: initializes StochasticForcing object with given parameters; 
/           de facto, this is a parametrized constructor;
/           since the parameters are not known prior to the declaration
/           of the object in Enzo, however, it is defined as a
/           regular method
/
/           Parameters:
/           my_spectral_rank -- dimension of Fouries space = grid rank
/           my_spect_profile -- shape of forcing power spectrum
/                               (1: delta peak, 2: band, 
/                                3: parabolic window) 
/           my_alpha -- ratio of domain length to integral length
/                       for each dimension (L = X/alpha)
/           my_band_width -- determines band width of the forcing 
/                            spectrum relative to alpha (maximal 
/                            value = 1)
/           my_intgr_vel -- characteristic velocity scale for each 
/                           dimension (charcteristic force per unit 
/                           mass F = V*V/L)
/           my_auto_corrl -- determines autocorrelation time of the
/                            stochastic force in units of the integral
/                            time scale T = L/V
/           my_soln_weight -- determines weight of solenoidal relative
/                             to dilatational modes (1 = purely 
/                             solenoidal, 0 = purely dilatational)
/
************************************************************************/


#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "StochasticForcing.h"

int StochasticForcing::Init(int my_spectral_rank,
			    forcing_type my_spect_profile, 
			    int *my_alpha, 
			    float *domain_length,
			    float *my_band_width,
			    float *my_intgr_vel,
		            float *my_auto_corrl,
		            float my_soln_weight,
                int my_seed)
{
    i1 = i2 = 0;
    j1 = j2 = 0;
    k1 = k2 = 0;
    decay = 0;

    AmpltThresh = 1.0e-4; // modes with normalised amplitude smaller than the threshold 
                          // are set equal to zero

    SpectralRank = my_spectral_rank;
    SpectProfile = my_spect_profile;

    /* set physical parameters */

    for (int dim = 0; dim < SpectralRank; dim++) {
	alpha[dim]         = my_alpha[dim];
	BandWidth[dim]     = (my_band_width[dim] > 0.0) ? my_band_width[dim] : 0.0;
	IntgrVelocity[dim] = my_intgr_vel[dim];
	IntgrLength[dim]   = (alpha[dim] > 0.0) ? (domain_length[dim] / alpha[dim]) : huge_number;
	IntgrTime[dim]     = (IntgrVelocity[dim] > 0.0) ? (IntgrLength[dim] / IntgrVelocity[dim]) : huge_number;
	AutoCorrlTime[dim] = my_auto_corrl[dim] * IntgrTime[dim];

	Amplitude[dim]     = NULL;
	InjectionEven[dim] = NULL;
	InjectionOdd[dim]  = NULL;
	SpectrumEven[dim]  = NULL;
	SpectrumOdd[dim]   = NULL;
    }


    /* determine boundary indices of the spectral domain */

    if (my_soln_weight >= 1.0) SolenoidalWeight = 1.0;
    if (my_soln_weight  < 1.0) SolenoidalWeight = my_soln_weight;
    if (my_soln_weight <= 0.0) SolenoidalWeight = 0.0;
    if (SpectralRank == 1) SolenoidalWeight = 0.0;
//TODO
/*
    if (MyProcessorNumber == ROOT_PROCESSOR) {
	cout << "Spectral profile " << SpectProfile << "\n"
	     << "Weight of solenoidal component " << SolenoidalWeight << "\n"
	     << "Alpha " << alpha[0] << " " << alpha[1] << " " << alpha[2] << " " << "\n"
	     << "Band width " << BandWidth[0] << " " << BandWidth[1] << " " << BandWidth[2] << " " << "\n"
	     << "Integral velocity " << IntgrVelocity[0] << " " << IntgrVelocity[1] << " " << IntgrVelocity[2] << " " << "\n"
	     << "Integral length " << IntgrLength[0] << " " << IntgrLength[1] << " " << IntgrLength[2] << " " << "\n"
	     << "Integral time " << IntgrTime[0] << " " << IntgrTime[1] << " " << IntgrTime[2] << " " << "\n"
	     << "Autocorrelation time " << AutoCorrlTime[0] << " " << AutoCorrlTime[1] << " " << AutoCorrlTime[2] << " " << "\n";
    }
*/
    /* determine boundary indices of spectral domain */

    if (SpectralRank > 0) {
	i1 = -2*alpha[0]+1; if (i1 > 0) i1 = 0;
	i2 =  2*alpha[0]-1; if (i2 < 0) i2 = 0;
    }
    if (SpectralRank > 1) {
	j1 = -2*alpha[1]+1; if (j1 > 0) j1 = 0;
	j2 =  2*alpha[1]-1; if (j2 < 0) j2 = 0;
    }
    if (SpectralRank > 2) {
	k1 = -2*alpha[2]+1; if (k1 > 0) k1 = 0;
	k2 =  2*alpha[2]-1; if (k2 < 0) k2 = 0;
    }

    /* determine number of linearly independent modes */
    
    NumModes = i2 + j2*(i2-i1+1) + k2*(j2-j1+1)*(i2-i1+1);
    
    mask = new int[NumModes];

    if (MyProcessorNumber == ROOT_PROCESSOR) {

	if (debug) printf("Total number of stochastic forcing modes = %"ISYM"\n",NumModes);

	/* determine forcing field amplitudes */

	for (int dim = 0; dim < SpectralRank; dim++) {
	    Amplitude[dim] = new float[NumModes];
	}

	if (SpectProfile == Peak) {

	    for (int n = 0; n < NumModes; n++)
		Amplitude[0][n] = 0.0;

	    Amplitude[0][(i2-1)/2] = 1.0;
	    Amplitude[0][(i2-1)/2+2] = 1.0;
	    if (SpectralRank > 1) 
		Amplitude[0][(j2+1)*(i2-i1+1)/2-1] = 1.0;
	    if (SpectralRank > 2) 
	    Amplitude[0][(k2+1)*(j2-j1+1)*(i2-i1+1)/2-1] = 1.0;

	} else {
	    
	    int i, j, k, n;
	    float x1, x2;
	    float a, a1 = 1.0, a2 = 1.0;
 
	    a1 = 1.0-BandWidth[0];
	    a2 = 1.0+BandWidth[0];
	    if (debug) {
		    printf("i1 = %"ISYM", i2 = %"ISYM"\n",i1,i2);
		    printf("a1 = %"FSYM", a2 = %"FSYM"\n\n",a1,a2);
	    }

	    /* compute amplitude factors for wave numbers within the interval [a1, a2] */

	    for (i = 1; i <= i2; i++) {
		a = 1.0; a = float(i) / float(alpha[0]);
		x1 = (a - a1); 
		x2 = (a2 - a);
		if (SpectProfile ==  Parabolic) {
		    Amplitude[0][i-1] = x1*x2;
		    if (Amplitude[0][i-1] <= 0.0) Amplitude[0][i-1] = 0.0;
		    Amplitude[0][i-1] *= Amplitude[0][i-1];
		} else if (SpectProfile == Band) {
		    if ((x1 >= 0.0) && (x2 >= 0.0)) Amplitude[0][i-1] = 1.0;
		}
		if (debug) {
		    printf("i = %"ISYM", a = %"FSYM"\n",i,a);
		    printf("x1 = %"FSYM", x2 = %"FSYM"\n",x1,x2);
		    printf("n = %"ISYM", Amplitude = %"FSYM"\n",i-1,Amplitude[0][i-1]);
		}
	    }

	    if (SpectralRank > 1) {
		
		float b, b1 = 1.0, b2 = 1.0;
		float f1 = 0.0, f2 = 0.0;
		float g1 = 1.0, g2 = 1.0;

		b1 = (a1 > 0.0) ? (1.0-BandWidth[1]) : 0.0; 
		b2 = (a1 > 0.0) ? (1.0+BandWidth[1]) : 2.0;
		if (b1 > 0.0) {
		    f1 = a1;
		    g1 = a1/b1; g1 *= g1;
		}
		if (b2 > 0.0) {
		    f2 = a2;
		    g2 = a2/b2; g2 *= g2;
		}
        /*
		if (debug) {
		    cout << "\nj1 = " << j1 << ", j2 = " << j2 << "\n";
		    cout << "b1 = " << b1 << ", b2 = " << b2 << "\n";
		    cout << "f1 = " << f1 << ", f2 = " << f2 << "\n";
		    cout << "g1 = " << g1 << ", g2 = " << g2 << "\n\n";
		}
*/
		/* compute amplitude factors for wave numbers bounded by the
		   ellipses with semi axes a1, b1 and a2, b2, respectively */
		
		n = i2;
		for (j = 1; j <= j2; j++) {
		    b = 0.0; if (alpha[1] > 0.0) b = float(j) / float(alpha[1]); b *= b;
		    for (i = i1; i <= i2; i++) {
			a = 0.0; if (alpha[0] > 0.0) a = float(i) / float(alpha[0]); a *= a;
			x1 = sqrt(a + g1*b) - f1;
			x2 = f2 - sqrt(a + g2*b);
			if (SpectProfile ==  Parabolic) {
			    Amplitude[0][n] = x1*x2;
			    if (Amplitude[0][n] <= 0.0) Amplitude[0][n] = 0.0;
			    Amplitude[0][n] *= Amplitude[0][n];
			} else if (SpectProfile == Band) {
			    if ((x1 >= 0.0) && (x2 >= 0.0)) Amplitude[0][n] = 1.0;
			}
            /*
			if (debug) {
			    cout << "i = " << i << ", a = " << a  << ", j = " << j << ", b = " << b << "\n";
			    cout << "x1 = " << x1 << ", x2 = " << x2 << "\n";
			    cout << "n = " << n << ", Amplitude = " << Amplitude[0][n] << "\n";
			}
            */
			++n;
		    }
		}

		if (SpectralRank > 2) {

		    float c, c1 = 1.0, c2 = 1.0;
		    float h1 = 1.0, h2 = 1.0;
		    
		    c1 = (b1 > 0.0) ? (1.0-BandWidth[2]) : 0.0;
		    c2 = (b1 > 0.0) ? (1.0+BandWidth[2]) : 0.0;
		    if (c1 > 0.0) {
			h1 = a1/c1; h1 *= h1;
		    }
		    if (c2 > 0.0) {
			h2 = a2/c2; h2 *= h2;
		    }
            /*
		    if (debug) {
			cout << "\nc1 = " << c1 << ", c2 = " << c2 << "\n";
			cout << "h1 = " << h1 << ", h2 = " << h2 << "\n\n";
		    }
		    */
		    /* compute amplitude factors for wave numbers bounded by the
		       ellipsoids with semi axes a1, b1, c1 and a2, b2, c2, respectively */

		    for (k = 1; k <= k2; k++) {
			c = 0.0; if (alpha[2] > 0.0) c = float(k) / float(alpha[2]); c *= c;
			for (j = j1; j <= j2; j++) {
			    b = 0.0; if (alpha[1] > 0.0) b = float(j) / float(alpha[1]); b *= b;
			    for (i = i1; i <= i2; i++) {
				a = 0.0; if (alpha[0] > 0.0) a = float(i) / float(alpha[0]); a*= a;
				x1 = sqrt(a + g1*b + h1*c) - f1;
				x2 = f2 - sqrt(a + g2*b + h2*c);
				if (SpectProfile ==  Parabolic) {
				    Amplitude[0][n] = x1*x2;
				    if (Amplitude[0][n] < 0.0) Amplitude[0][n] = 0.0;
				    Amplitude[0][n] *= Amplitude[0][n];
				} else if (SpectProfile == Band) {
				    if ((x1 >= 0.0) && (x2 >= 0.0)) Amplitude[0][n] = 1.0;
				}
                /*
				if (debug) {
				    cout << "i = " << i << ", a = " << a  << ", j = " << j << ", b = " << b 
					 << ", k = " << k << ", c = " << c << "\n";
				    cout << "x1 = " << x1 << ", x2 = " << x2 << "\n";
				    cout << "n = " << n << ", Amplitude = " << Amplitude[0][n] << "\n";
				}
                */
				++n;
			    }
			}
		    }
		}
	    }
	}

	/* normalise amplitude factors and 
           set flags for modes with amplitude larger than the threshold */

	float norm = 0.0;
	for (int n = 0; n < NumModes; n++) norm += 2*Amplitude[0][n]*Amplitude[0][n];
	norm = 1/sqrt(norm);

	NumNonZeroModes = 0;
	for (int n = 0; n < NumModes; n++) {
	    Amplitude[0][n] *= norm;
	    if (Amplitude[0][n] > AmpltThresh) {
		mask[n] = 1;
		++NumNonZeroModes;
	    } else {
		mask[n] = 0;
	    }
	    if (SpectralRank > 1) Amplitude[1][n] = Amplitude[0][n];
	    if (SpectralRank > 2) Amplitude[2][n] = Amplitude[1][n];
	    //if (debug) cout << n << "   " << Amplitude[0][n] << " " << mask[n] << "\n";
	}	    	    
	printf("Number of non-zero stochastic forcing modes = %"ISYM"\n",NumNonZeroModes);
	
	if (NumNonZeroModes == 0) return FAIL;
	
	for (int dim = 0; dim < SpectralRank; dim++) {

	    InjectionEven[dim] = new float[NumModes];
	    InjectionOdd[dim]  = new float[NumModes];

	    norm = 3.0 * IntgrVelocity[dim] * IntgrVelocity[dim] / 
		 (sqrt(1.0 - 2.0*SolenoidalWeight + 3.0*SolenoidalWeight*SolenoidalWeight) * IntgrLength[dim]);

//	    if (debug) cout << "Normalization factor[" << dim << "] = " << norm << "\n";

	    for (int n = 0; n < NumModes; n++) {
		Amplitude[dim][n] *= norm;
		InjectionEven[dim][n] = 0.0;
		InjectionOdd [dim][n] = 0.0;
	    }
        }

	/* initialise new sequence of random numbers */

	seed = my_seed;
//	if (debug) cout << "Seed of random number generator = " << seed << "\n";
	RandUni(seed);

	/* compute initial set of random deviates */

	Inject();
    }

    /* communicate mask among processors */

    CommunicationBroadcastFlags();

    /* determine number of non-zero modes if not root */

    if (MyProcessorNumber != ROOT_PROCESSOR) {
	NumNonZeroModes = 0;
	for (int n = 0; n < NumModes; n++)
	    if (mask[n]) ++NumNonZeroModes;
//	if (debug) cout << "Number of non-zero stochastic forcing modes = " << NumNonZeroModes 
//			<< ", proc #" << MyProcessorNumber << "\n";
    }

    /* allocate memory for the forcing spectrum and set inital random deviates */

    for (int dim = 0; dim < SpectralRank; dim++) {

	SpectrumEven[dim] = new float[NumNonZeroModes];
	SpectrumOdd[dim]  = new float[NumNonZeroModes];
	
	if (MyProcessorNumber == ROOT_PROCESSOR) {
	    for (int n = 0, m = 0; n < NumModes; n++) {
		if (mask[n]) {
		    SpectrumEven[dim][m] = InjectionEven[dim][n];
		    SpectrumOdd [dim][m] = InjectionOdd [dim][n];
		    ++m;
		}
	    }
	}
    }

    /* communicate spectrum among processors */

    CommunicationBroadcastSpectrum();

    return SUCCESS;
}
