
/**************************************************************************
 *
 *  STOCHASTIC FORCING CLASS
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1:
 *
 *  PURPOSE: composes and evolves a stochastic driving force in Fourier space
 *
 **************************************************************************/
#ifndef STOCHASTIC_FORCING_DEFINED__
#define STOCHASTIC_FORCING_DEFINED__
//#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"

class StochasticForcing {

 private:

    int SpectralRank;                   // rank of the spectrum
    int NumModes;                       // number of Fourier modes in the spectral domain
    int NumNonZeroModes;                // number of non-zero Fourier modes
    int i1, i2, j1, j2, k1, k2;         // index boundaries specifying the range of wave numbers
    int decay;                          // if set to non-zero value, the force field decays
    float AmpltThresh;                  // threshold value for the normalised amplitude of modes counted as non-zero
//
//  Integral scale parameters (these are generally NOT the expectation values but characteristic quantities!)
//
    int alpha[MAX_DIMENSION];           // ratio of domain size to integral scale
    float BandWidth[MAX_DIMENSION];     // bandwidth of the forcing spectrum in units of 2*alpha
    float IntgrVelocity[MAX_DIMENSION]; // integral velocity V in the statistically stationary regime
    float IntgrTime[MAX_DIMENSION];     // integral (large-eddy turn-over) time T=L/V of the flow evolution
    float IntgrLength[MAX_DIMENSION];   // integral length L of the force field
    float WaveNumber[MAX_DIMENSION];    // characteristic wave number = 2 pi/L
//  
//  Stochastic process parameters (in Fourier space)
// 
    forcing_type SpectProfile;          // profile of forcing spectrum
    float AutoCorrlTime[MAX_DIMENSION]; // autocorrelation time of the stochastic process
    float SolenoidalWeight;             // determines weight of solenoidal relative to dilatational components
    float DecayInitTime;                // time at which the decay of the force field is initiated
//
//  Spectral data
//
    int *mask;                           // flags identifying non-zero modes
    float *Amplitude[MAX_DIMENSION];     // amplitudes of forcing modes
    float *InjectionEven[MAX_DIMENSION]; // random increments (cos modes)
    float *InjectionOdd[MAX_DIMENSION];  // random increments (sin modes)
    float *SpectrumEven[MAX_DIMENSION];  // forcing cos modes
    float *SpectrumOdd[MAX_DIMENSION];   // forcing sin modes

 public:

//
// Constructor (nullifies data)
//
    StochasticForcing();
//
// Set forcing amplitudes and initialise spectrum
//
    int Init(int my_spectral_rank,
	     forcing_type my_spect_profile, 
	     int *my_alpha, 
	     float *domain_length,
	     float *my_band_width_fct,
	     float *my_intgr_vel,
	     float *my_auto_corrl_time,
	     float my_soln_weight,
       int my_seed);
//
// Destructor
//
    ~StochasticForcing();
//
// Evolve the force field over a time step which is small compared to the autocorrelation time
//
    void Evolve(float dt);
//
// Calculate the instantaneous RMS force
//
    float RMS(void);
//
// Write/read forcing spectrum to/from output file
//
    int ReadSpectrum(char *fname);
    int WriteSpectrum(char *fname);
    void WriteParameters(FILE *fptr);
// 
// Initiate the decay of the force field
//
    void set_decay(void);
//
// Set the weighting parameter
//
    void set_SolenoidalWeight(int my_soln_weight);
//
// Get boundary indices and range of spectral domain
//
    int get_LeftBoundary(int dim);
    int get_RightBoundary(int dim);
    int get_Range(int dim);
//
// Get number of forcing modes
//
    int get_NumModes(void);
    int get_NumNonZeroModes(void);
//
// Get spectral profile
//
    forcing_type get_SpectProfile(void);
//
// Get the characteristic wave number
//
    float get_WaveNumber(int dim);
//
// Get the integral length scale
//
    float get_IntgrLength(int dim);
//
// Get the integral time scale
//
    float get_IntgrTime(int dim);
//
// Copy flags identifying non-zero modes
//
    void copy_mask(int* target);
//
// Copy non-zero modes of the forcing spectrum
//
    void copy_SpectrumEven(int dim, float* target);
    void copy_SpectrumOdd(int dim, float* target);
//
// Copy complete forcing spectrum into a single array for subsequent processing
//
    void copy_ExpandedSpectrum(int dim, float* target);

 private:

//
// Compute random injection
//
    void Inject(void);
//
// Random number generator
//
    void GaussDeviate(float amplt, float *x, float *y);
    double RandUni(int &idum);

//
// Communicate forcing spectrum from root to other processors
//
    void CommunicationBroadcastSpectrum(void);
//
// Communicate flags identifying non-zero modes from root to other processors
//
    void CommunicationBroadcastFlags(void);
};


inline void StochasticForcing::set_decay(void) 
{
    decay = 1;
}

inline void StochasticForcing::set_SolenoidalWeight(int my_soln_weight)
{
    if (my_soln_weight >= 1.0) SolenoidalWeight = 1.0;
    if (my_soln_weight  < 1.0) SolenoidalWeight = my_soln_weight;
    if (my_soln_weight <= 0.0) SolenoidalWeight = 0.0;
}

inline int StochasticForcing::get_LeftBoundary(int dim)
{
    switch (dim) {
	case 1:
	    return i1;
	case 2:
	    return j1;
	case 3:
	    return k1;
	default:
	    return 0;
    }
}

inline int StochasticForcing::get_RightBoundary(int dim)
{
    switch (dim) {
	case 1:
	    return i2;
	case 2:
	    return j2;
	case 3:
	    return k2;
	default:
	    return 0;
    }
}

inline int StochasticForcing::get_NumModes(void)
{
    return NumModes;
}

inline int StochasticForcing::get_NumNonZeroModes(void)
{
    return NumNonZeroModes;
}

inline forcing_type StochasticForcing::get_SpectProfile(void)
{
    return SpectProfile;
}

inline float StochasticForcing::get_WaveNumber(int dim)
{
    return WaveNumber[dim];
}

inline float StochasticForcing::get_IntgrLength(int dim)
{
    return IntgrLength[dim];
}

inline float StochasticForcing::get_IntgrTime(int dim)
{
    return IntgrTime[dim];
}

inline void StochasticForcing::copy_mask(int* target)
{
    for (int n = 0; n < NumModes; n++) 
	target[n] = mask[n];
}

inline void StochasticForcing::copy_SpectrumOdd(int dim, float* target)
{
    for (int m = 0; m < NumNonZeroModes; m++) 
	target[m] = SpectrumOdd[dim][m];
}

inline void StochasticForcing::copy_SpectrumEven(int dim, float* target)
{
    for (int m = 0; m < NumNonZeroModes; m++) 
	target[m] = SpectrumEven[dim][m];
}

inline void StochasticForcing::copy_ExpandedSpectrum(int dim, float* target)
{
    for (int n = 0, m = NumNonZeroModes-1; n < NumModes; n++) 
	target[n] = mask[NumModes-1-n] ? SpectrumOdd[dim][m--] : 0.0;

    target[NumModes] = 0.0;

    for (int n = NumModes+1, m = 0; n <= 2*NumModes; n++) 
	target[n] = mask[n-NumModes-1] ? SpectrumEven[dim][m++] : 0.0;
}

#endif
