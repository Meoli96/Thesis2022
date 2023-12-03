#include "../include/pch.h"
#include "limits.h"

#include "../include/HarmonicAcc.h"
#include "../include/moon.h"
#include "../include/utils.h"

std::string dataPath = "data/";
std::string cofDir = dataPath + "cof/";
std::string radDir = dataPath + "rad/";

int main(int argc, char *argv[])
{

    // Argv[1] specifies the degree of the harmonic expansion
    int degree = 0;
    std::string str_deg;
    size_t N = 0;
    std::chrono::microseconds duration{};

    if (argc == 2)
    {

        degree = std::stoi(argv[1]);
        str_deg = std::string{argv[1]};
    }

    if (degree != 0 && !str_deg.empty())
    {
        try
        {
            
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        
        // Create harmonic filepath and radius path
        std::string cofFile = cofDir + "cilm" + str_deg + ".csv";
        std::string radFile = radDir + "r360.txt";

        std::ofstream f{"acc.csv"};
        // Initialize HarmonicAcc class
        HarmonicAcc moon{degree, cofFile};

        // Now load radius from file
        double *r;
        loadRadiusGMAT(r, radFile, N);
        double *a = new double[3 * N];

        // Now I can compute the acceleration
        for (size_t i = 0; i < N; i++)
        {
            // Get r pointer
            double *curr_r = &r[i * 3];
            double *curr_a = &a[i * 3];

            // Compute point mass acceleration
            double a_dsh[3];
            double a_pm[3];

            auto start = std::chrono::high_resolution_clock::now();
            moon.getPointAcc(curr_r, a_pm, muun);
            moon.getAcc(curr_r, a_dsh, muun, r_cmoon);
            auto stop = std::chrono::high_resolution_clock::now();
            curr_a[0] = a_pm[0] + a_dsh[0];
            curr_a[1] = a_pm[1] + a_dsh[1];
            curr_a[2] = a_pm[2] + a_dsh[2];

            duration += std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            f << curr_a[0] << "," << curr_a[1] << "," << curr_a[2] << std::endl;
        }
   

    
    }

    else
    {
        std::cout << "No degree provided. Exiting" << std::endl;
    }
}