// Compute disturbing acceleration usign Pines form
#include "../include/HarmonicAcc.h"
#include "../include/utils.h"

//////////////////////////////////////// CONSTRUCTOR /////////////////////////////////////////

HarmonicAcc::HarmonicAcc(uint32_t degree, std::string path)
{
    /*
        Constructor for HarmonicAcc class.
        Is responsible for allocating and managing memory needed
        to store all whats needed for disturbing acceleration computation.

        __NVCC__ VERSION: This constructor should allocate on the target GPU
            - cilm -> also to copy on construction
            - Alm
            - reim
        and store the device pointers of each memory segment.
    */

    // Alloc for cilm and get them from path

    this->path = path;
    this->degree = degree;
    this->pow2 = ceil2pow(degree);

    // Compute this->bufsize and padding
    size_t L = this->degree * (this->degree + 1) / 2 + this->degree;
    this->padding = (32 - L % 32);
    this->bufsize = L + padding;

    this->cilm = new double[2 * this->bufsize]();
    this->Alm = new double[this->bufsize + degree + 2]();
    this->reim = new double[2 * (degree + 1)]();

    // Allocate for legendre functions
    /*
        Alm needs to be a [degree+1, degree+1] lower triangular matrix,
        so im going to allocate a [degree, degree] with this->bufsize and
        then add another row of 361 terms with [degree + 2]

        This is because  D(degree, degree) = ...*Alm(degree+1,degree+1)
    */

    this->N1 = new double[this->bufsize + degree + 2]();
    this->N2 = new double[this->bufsize + degree + 2]();

    // Alloc normalization correction matrices
    this->V01 = new double[bufsize]();
    this->V11 = new double[bufsize]();

    // Alloc sums buffer
    this->sums = new double[4];

    // Read harmonic coefficients from path
    readCSV(this->path, this->cilm, 2 * this->bufsize, this->padding);

    // Precompute diagonal elements of Alm (not dependent on u)
    Alm[0] = 1.0;

    for (uint16_t l = 1; l <= degree + 1; ++l)
    {

        size_t idx00 = ij2idx(l, l);
        size_t idx11 = idx00 - (l + 1);

        // A(n,n) = sqrt( ( 2*n + 1 ) / (2 * n) ) * A(n-1, n-1)
        Alm[idx00] = sqrt(double(2 * l + 1) / double(2 * l)) * Alm[idx11];
    }

    // Compute N1 and N2

    for (uint16_t m = 0; m <= degree + 1; ++m)
    {
        for (uint16_t l = m + 2; l <= degree + 1; ++l)
        {
            size_t idx = ij2idx(l, m);
            N1[idx] = sqrt(double((2 * l + 1) * (2 * l - 1)) / double((l + m) * (l - m)));
            N2[idx] = sqrt(double((2 * l + 1) * (l - m - 1) * (l + m - 1)) /
                           double((2 * l - 3) * (l + m) * (l - m)));
        }
    }

    // Compute V01, V11
    for (int l = 0; l <= degree; l++)
    {
        for (int m = 0; m <= l; m++)
        {
            size_t idx = ij2idx(l, m);
            V01[idx] = sqrt(double((l - m) * (l + m + 1)));
            V11[idx] = sqrt(double((2 * l + 1) * (l + m + 2) * (l + m + 1)) / double((2 * l + 3)));
            if (m == 0)
            {
                V01[idx] /= sqrt(double(2));
                V11[idx] /= sqrt(double(2));
            }
        }
    }
}

/////////////////////////////////////////// DESTRUCTOR ///////////////////////////////////////

HarmonicAcc::~HarmonicAcc()
{

    delete[] this->sums;

    delete[] this->cilm;
    delete[] this->Alm;
    delete[] this->reim;
}

////////////////////////////////////////////// __NVCC__ KERNEL //////////////////////////////////////

/////////////////////////////////////// getAcc ////////////////////////////////////////////

void HarmonicAcc::getAcc(double *r, double *acc, double mi, double r_cb)
{
    /// @brief
    /// @param r
    //      Vector radius [r.x, r.y, r.z]
    /// @param mi
    //      Gravitational parameter of central body
    /// @param a
    //      Central body equatorial radius (or mean radius)
    /// @return acc
    //

    double rmod = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    s = r[0] / rmod;
    t = r[1] / rmod;
    u = r[2] / rmod;

    // Initialize sum terms
    for (int i = 0; i < 4; i++)
    {
        this->sums[i] = 0;
    }

    double rho = r_cb / rmod;
    double rho1 = -(mi / rmod) * rho;

    // Initialize reim
    this->reim[0] = 1;
    this->reim[1] = 0;
    // Initialize rho_n

    // Compute reim
    for (int i = 1; i <= this->degree; i++)
    {
        // r_m = s*r_m-1 - t*i_m-1
        this->reim[i * 2] =
            s * this->reim[i * 2 - 2] - t * this->reim[i * 2 - 1];
        // i_m = s*i_m-1 + t*r_m-1
        this->reim[i * 2 + 1] =
            s * this->reim[i * 2 - 1] + t * this->reim[i * 2 - 2];
    }
    // Compute ALF

    // A(1,0)
    Alm[1] = u * sqrt(3.0);

    for (int l = 1; l <= degree; l++)
    {
        // Off diagonal elements
        size_t idx00 = ij2idx(l, l);  // (l, m) -> idx
        size_t idx10 = idx00 + l + 1; // (l+1, m) -> idx

        Alm[idx10] = u * sqrt(double(2 * l + 3)) * Alm[idx00];
    }
    for (int m = 0; m <= degree + 1; m++)
    {
        for (int l = m + 2; l <= degree + 1; l++)
        {
            // Column recursion algorithm
            size_t idx00 = ij2idx(l, m);
            size_t idx10 = idx00 - l;     //(l-1, m)
            size_t idx20 = idx10 - l + 1; //(l-2, m)
            Alm[idx00] = u * N1[idx00] * Alm[idx10] -
                         N2[idx00] * Alm[idx20];
        }
    }


    double sqrt2 = sqrt(2);
    for (uint16_t l = 1; l <= degree; l++)
    {
        // Compute rho1
        rho1 *= rho;
        // Initialize sum terms in m
        double temp0 = 0;
        double temp1 = 0;
        double temp2 = 0;
        double temp3 = 0;
        for (uint16_t m = 0; m <= l; m++)
        {
            // Compute linear index
            size_t idx00 = ij2idx(l, m);  // (l,m) -> idx
            size_t idx11 = idx00 + l + 2; // (l+1,m+1) -> idx
            size_t idx01 = idx00 + 1;     // (l,m+1 -> idx)

            double C = cilm[2 * idx00];
            double S = cilm[2 * idx00 + 1];

            // D(l, m) = C(l,m)*r(m) + S(l,m)*i(m)
            double Dlm = (C * reim[2 * m] + S * reim[2 * m + 1]) * sqrt2;

            /*

                Here I need to check if m==0 (with ternary operator)
                or else reim goes out of bounds :(

            */

            // E(l, m) = C(l,m)*r(m-1) +  S(l, m)*i(m-1)
            double Elm = m == 0 ? 0 : (C * reim[2 * m - 2] + S * reim[2 * m - 1]) * sqrt2;

            // F(l, m) = S(l,m)*r(m-1) - C(l,m)*i(m-1)
            double Flm = m == 0 ? 0 : (S * reim[2 * m - 2] - C * reim[2 * m - 1]) * sqrt2;

            // Correct normalization for Legendre coefficients
            /*
                If l==m A(l,l+1) = d (A(l,l)/du) == 0
            */
            double A00 = Alm[idx00];
            double A01 = m == l ? 0 : Alm[idx01] * V01[idx00];
            double A11 = Alm[idx11] * V11[idx00];

            // Compute partial sum

            // temp0 = sum_m Elm(l,m)*A(l,m), s component
            temp0 += m * Elm * A00;
            // sum1 = sum_m Flm(l,m)*A(l,m), t component
            temp1 += m * Flm * A00;
            // temp2 = sum_m D(l,m)*A(l,m+1), u component
            temp2 += Dlm * A01;
            // temp3 = sum_m D(l,m)*A(l+1, m+1)*c2, R component
            temp3 += Dlm * A11;
        }
        // Multiply by rho
        double fact = rho1 / r_cb;
        this->sums[0] += fact * temp0;
        this->sums[1] += fact * temp1;
        this->sums[2] += fact * temp2;
        this->sums[3] -= fact * temp3;
    }
    // Get acceleration

    // i-component
    acc[0] = this->sums[0] + this->sums[3] * s;
    // j-component
    acc[1] = this->sums[1] + this->sums[3] * t;
    // k-component
    acc[2] = this->sums[2] + this->sums[3] * u;

    return;
}

///////////////////////////////////////// getPointAcc /////////////////////////////////////

void HarmonicAcc::getPointAcc(double *r, double *acc, double mi)
{
    double rmod = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double fact = -mi / (rmod * rmod * rmod);
    acc[0] = fact * r[0];
    acc[1] = fact * r[1];
    acc[2] = fact * r[2];
}
