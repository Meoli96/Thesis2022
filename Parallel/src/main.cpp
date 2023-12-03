#include "pch.h"
#include "limits.h"
#include "fstream"
#include "utils.h"

#include <chrono>

// Global variable declarations (fixed constants of the problem)
constexpr double muun = 4904.8695; // km^3/s^2
constexpr double r_cmoon = 1737.4; // km, moon radius

struct cuArgs
{
    int degree;
    size_t bufsize;
    size_t pow2;
    double sqrt2;
};

__global__ void getAdsh(double *cilm, double *Alm, double *reim, double *parts, struct cuArgs *d_cu)
{
    /// @brief This function computes the disturbing acceleration on the GPU
    /// @param cilm Pointer to Clm and Slm formatted as [ C(l,m), S(l,m), C(l,m+1), S(l,m+1) ... ]
    /// @param Alm Pointer to the legendre matrix
    /// @param reim Pointer to the Pines variables formatted as [ re(m), im(m), re(m+1), im(m+1), ... ]
    /// @param parts Pointer to the partial sums buffer
    /// @param bufsize multiple of size of various buffers
    /// @param pow2 Closer power of 2 to the degree of the problem
    /// @return
    extern __shared__ double shMem[];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int shid = tid;

    int mul;

    // Get additional arguments
    int degree = d_cu->degree;
    size_t bufsize = d_cu->bufsize;
    size_t pow2 = d_cu->pow2;
    double sqrt2 = d_cu->sqrt2;

    size_t buf_stop = 4 * bufsize;

    // Need to spawn at least 4*pow2 threads to initialize shared memory
    if (tid < 4 * pow2)
    {
        shMem[tid] = 0;
    }
    __syncthreads();

    while (tid < 4 * bufsize)
    {

        if (tid < bufsize)
        {
            mul = 0;
        }
        if (tid >= bufsize && tid < bufsize * 2)
        {
            mul = 1;
        }
        if (tid >= bufsize * 2 && tid < bufsize * 3)
        {
            mul = 2;
        }
        if (tid >= bufsize * 3 && tid < bufsize * 4)
        {
            mul = 3;
        }

        // idx range is [0, bufsize]
        int idx = tid - mul * bufsize;

        int l = floor(-0.5 + sqrt(0.25 + 2 * idx));
        int m = idx - l * (l + 1) / 2;

        // Initialize shared memory to 0

        double C = cilm[2 * idx];
        double S = cilm[2 * idx + 1];

        /*
            I need to account for padding bytes. Those will produce wrong values of l
            causing out of bounds memory access
        */
        if (l <= degree)
        {

            if (mul == 0)
            {
                /*
                   This if statement shouldnt (and i capitalize it because i really dont know if thats the case)
                   cause warp divergency, because bufsize % 32 = 0 and thus multiple of warpsize
               */

                // Access A00
                double A00 = Alm[idx];
                // Compute Elm

                double re2 = m == 0 ? 0 : reim[2 * m - 2];
                double im2 = m == 0 ? 0 : reim[2 * m - 1];

                // E(l, m) = C(l,m)*r(m-1) +  S(l, m)*i(m-1)

                double Elm = (C * re2 + S * im2) * sqrt2;
                double fact = m * Elm * A00;
                atomicAdd(&shMem[l + mul * pow2], fact);
            }
            else if (mul == 1)
            {

                double A00 = Alm[idx];

                double re2 = m == 0 ? 0 : reim[2 * m - 2];
                double im2 = m == 0 ? 0 : reim[2 * m - 1];

                // F(l, m) = S(l,m)*r(m-1) - C(l,m)*i(m-1)
                double Flm = (S * re2 - C * im2) * sqrt2;
                double fact = m * Flm * A00;
                atomicAdd(&shMem[l + mul * pow2], fact);
            }
            else if (mul == 2)
            {
                // Compute Dlm
                size_t idx01 = idx + 1;
                // Access A01
                double A01 = m == l ? 0 : Alm[idx01];
                // Compute Dlm*Plm(l, m+1);
                double re = reim[2 * m];
                double im = reim[2 * m + 1];

                double c1 = sqrt(double((l - m) * (l + m + 1)));
                if (m == 0)
                {
                    c1 /= sqrt2;
                }

                double Dlm = (C * re + S * im) * sqrt2;
                double fact = Dlm * A01 * c1;
                atomicAdd(&shMem[l + mul * pow2], fact);
            }
            else if (mul == 3)
            {
                size_t idx11 = idx + l + 2;
                // Get Legendre value
                double A11 = Alm[idx11];

                // Get values from memory
                double re = reim[2 * m];
                double im = reim[2 * m + 1];

                // Compute Dlm*Plm(l+1, m+1)

                double c2 = sqrt(double((l + m + 2) * (l + m + 1)) /
                                 double((2 * l + 3) * (2 * l + 2)));

                if (m == 0)
                {
                    c2 /= sqrt2;
                }

                // D(l, m) = C(l,m)*r(m) + S(l,m)*i(m)
                double Dlm = (C * re + S * im) * sqrt2;
                double fact = Dlm * A11 * c2;
                atomicAdd(&shMem[l + mul * pow2], fact);
            }
        }
        // Syncthreads here, load shared memory back to global memory

        tid += gridDim.x;
    }
    __syncthreads();
    if (shid < 4 * pow2)
    {
        parts[shid] += shMem[shid];
    }
}

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

    HarmonicAcc(uint32_t degree, std::string path)
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

        cudaError_t cudaErr;

        this->path = path;
        this->degree = degree;
        this->pow2 = ceil2pow(degree);

        // Compute this->bufsize and padding
        size_t L = this->degree * (this->degree + 1) / 2 + this->degree;
        this->padding = (32 - L % 32);
        this->bufsize = L + padding;
        // Alloc cuArgs structure
        struct cuArgs cu
        {
            degree, bufsize, pow2, sqrt(double(2))
        };
        /*
            Each buffer that is interacting with the GPU will be allocated
            as page-locked host memory to increase performance of transfer

        */
        cudaErr = cudaMallocHost(&(this->Alm), sizeof(double) * (bufsize + degree + 2));
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMallocHost(&(this->cilm), sizeof(double) * bufsize * 2);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMallocHost(&(this->reim), sizeof(double) * 2 * (degree + 1));
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

        // Host result buffer
        cudaErr = cudaMallocHost(&(this->h_sums), sizeof(double) * pow2 * 4);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

        // Allocate Alm, cilm and reim, cuArgs on device
        cudaErr = cudaMalloc((void **)(&d_Alm), sizeof(double) * (bufsize + degree + 2));
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMalloc((void **)(&d_reim), sizeof(double) * 2 * (degree + 1));

        cudaErr = cudaMalloc((void **)(&d_cilm), sizeof(double) * (bufsize * 2));
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMalloc((void **)(&d_sums), sizeof(double) * pow2 * 4);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMemset(d_sums, 0, sizeof(double) * pow2 * 4);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMalloc((void **)&d_cu, sizeof(struct cuArgs));
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

        // Allocate for legendre functions
        /*
            Alm needs to be a [degree+1, degree+1] lower triangular matrix,
            so im going to allocate a [degree, degree] with this->bufsize and
            then add another row of 361 terms with [degree + 2]

            This is because  D(degree, degree) = ...*Alm(degree+1,degree+1)
        */

        this->N1 = new double[this->bufsize + degree + 2]();
        this->N2 = new double[this->bufsize + degree + 2]();

        // Alloc sums buffer
        this->sums = new double[4];

        // Read harmonic coefficients from path
        readCSV(this->path, this->cilm, 2 * this->bufsize, this->padding);

        // Memcpy harmonic coefficients to device memory
        cudaErr = cudaMemcpy(d_cilm, cilm, sizeof(double) * 2 * bufsize, cudaMemcpyHostToDevice);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }
        cudaErr = cudaMemcpy(d_cu, &cu, sizeof(struct cuArgs), cudaMemcpyHostToDevice);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

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
    }

    // Destructor
    ~HarmonicAcc()
    {

        delete[] this->sums;

        // Deallocate host memory
        cudaFreeHost(this->cilm);
        cudaFreeHost(this->Alm);
        cudaFreeHost(this->reim);
        cudaFreeHost(this->h_sums);
    }

    // Getters
    size_t getBufsize() { return this->bufsize; }

    void getAcc(double *r, double *acc, double mi, double r_cb)
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

        // __NVCC__ code

        cudaError_t cudaErr;

        // memcpy reim and Pilm

        cudaErr = cudaMemcpy(d_reim, reim, sizeof(double) * 2 * (degree + 1), cudaMemcpyHostToDevice);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

        cudaErr = cudaMemcpy(d_Alm, Alm, sizeof(double) * (bufsize + degree + 2), cudaMemcpyHostToDevice);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

        // Call kernel
        // Get number of threads per block and blocks per grid

        size_t shSize = sizeof(double) * pow2 * 4;
        getAdsh<<<1024, 512, shSize>>>(d_cilm, d_Alm, d_reim, d_sums, d_cu);
        cudaErr = cudaGetLastError();
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
        }
        cudaErr = cudaDeviceSynchronize();
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
        }
        // memcpy results back to host buf
        cudaErr = cudaMemcpy(h_sums, d_sums, sizeof(double) * pow2 * 4, cudaMemcpyDeviceToHost);
        if (cudaErr != cudaSuccess)
        {
            std::cout << cudaGetErrorString(cudaErr);
            throw;
        }

        double fact = 0;
        // Now iterate over h_sums and compute
        for (int l = 0; l <= degree; l++)
        {
            // For l==0 we dont update rho1, we also make sure
            // that l == 0 doesnt contribute to the sum
            if (l)
            {
                rho1 *= rho;
                fact = rho1 / r_cb;
            }

            // Elm
            this->sums[0] += fact * h_sums[l];
            // Flm
            this->sums[1] += fact * h_sums[l + pow2];
            // Dlm * A01
            this->sums[2] += fact * h_sums[l + 2 * pow2];
            // Dlm * A11
            this->sums[3] -= fact * h_sums[l + 3 * pow2];
        }

        // i-component
        // i-component
        acc[0] = this->sums[0] + this->sums[3] * s;
        // j-component
        acc[1] = this->sums[1] + this->sums[3] * t;
        // k-component
        acc[2] = this->sums[2] + this->sums[3] * u;

        return;
    }
    void getPointAcc(double *r, double *acc, double mi)
    {
        double rmod = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        double fact = -mi / (rmod * rmod * rmod);
        acc[0] = fact * r[0];
        acc[1] = fact * r[1];
        acc[2] = fact * r[2];

        return;
    }

private:
    std::string path; // Path to harmonic coefficients
    size_t degree;    // Degree of harmonic coefficients

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

    /* Device buffers */
    double *d_cilm;
    double *d_Alm;
    double *d_reim;
    struct cuArgs *d_cu;
    // Input buffers
    double *d_sums;
    double *h_sums;
};

std::string dataPath = "/home/peppe/Tesi/code/Parallel/data/";
std::string cofDir = dataPath + "cof/";
std::string radDir = dataPath + "rad/";

int main(int argc, char *argv[])
{
    std::chrono::microseconds duration{};
    // Argv[1] specifies the degree of the harmonic expansion
    uint16_t degree = 0;
    std::string str_deg;
    size_t N = 0;
    int temp = 0;
    std::ofstream f{"time.txt", std::ios::app};

    if (argc == 2)
    {
        degree = std::stoi(argv[1]);
        str_deg = std::string{argv[1]};
    }

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    if (degree != 0 && !str_deg.empty())
    {
        try
        {
            // Create harmonic filepath and radius path
            std::string cofFile = cofDir + "cilm" + str_deg + ".csv";
            std::string radFile = radDir + "r360.txt";

            // Initialize HarmonicAcc class
            HarmonicAcc moon{degree, cofFile};

            // Now load radius from file
            double *r;
            loadRadiusGMAT(r, radFile, N);

            // Allocate acceleration buffer
            double *a = new double[N * 3];

            // Now I can compute the acceleration
            for (size_t i = 0; i < 3; i++)
            {
                // Get r pointer
                double *curr_r = &r[i * 3];
                // Get a pointer
                double *curr_a = &a[i * 3];

                curr_a[0] = 0;
                curr_a[1] = 0;
                curr_a[2] = 0;

                // Compute point mass acceleration
                double a_dsh[3];

                auto start = std::chrono::high_resolution_clock::now();
                moon.getAcc(curr_r, a_dsh, muun, r_cmoon);
                auto stop = std::chrono::high_resolution_clock::now();

                curr_a[0] = a_dsh[0];
                curr_a[1] = a_dsh[1];
                curr_a[2] = a_dsh[2];
                temp = i;
                duration += std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            }
        }
        catch (...)
        {
            duration /= 3;
            f << duration.count() << ", " << degree << std::endl;
            return 0;
        }
        std::cout << duration.count() << ", " << degree << std::endl;
    }

    else
    {
        std::cout << "No degree provided. Exiting" << std::endl;
    }
}