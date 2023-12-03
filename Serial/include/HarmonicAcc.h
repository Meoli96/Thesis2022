#include "pch.h"

#ifdef __NVCC__
__global__ void getAdsh(double *, double *,
                        double *, double *, double *, double *);
#endif

class HarmonicAcc
{

    /*
        This singleton class will be responsible of holding all of the data and
       also the computing of the disturbing acceleration given by spherical
       harmonics.

        It will be responsible of the inizialization from a csv file of the Clm
       and Slm coefficients, and also it will manage all the memory allocation
       for the variables needed for the computing.

    */
public:
    // Constructors
    HarmonicAcc();

    HarmonicAcc(uint32_t, std::string);
    // Destructor
    ~HarmonicAcc();

    // Getters
    size_t getBufsize() { return this->bufsize; }

    void getAcc(double *r, double *acc, double mi, double a);
    void getPointAcc(double *, double *, double);

private:
    std::string path; // Path to harmonic coefficients
    uint32_t degree;  // Degree of harmonic coefficients

    // Buffer dimensions and padding
    size_t bufsize;
    char padding;

    // Various buffers
    double *cilm; // Harmonic coefficients

    // Stuff needed for associated legendre functions
    double *Alm; // ALF
    // Normalizing matrices
    double *N1;
    double *N2;

    // Pines variable
    

    double *reim; // Real and Imaginary part for expansion
    size_t pow2;  // Upper closer power of 2 to degree
 
    double s, t, u;

    // Array holding all temporary summatories
    double *sums;

#ifdef __NVCC__
    /* Device buffers */
    // Input buffers
    double *d_cilm;
    double *d_Alm;
    double *d_reim;
    // Output buffer
    double *d_sums;
    double *h_sums;

#else

    // Normalization correction
    /*
        Used to correct normalization of A(l, m+1) and
        A(l+1, m+1) respectively.
    */
    double *V01;
    double *V11;

#endif
};