#include "../include/utils.h"

void readCSV(std::string path, double *buf, size_t bufsize, uint8_t padding)
{

    /*
        Clm and Slm are going to be stored contigously, eg [C11, S11, C20, S20,
       ...]

        This gives advantages in both
        serial and parallel case.

        Legendre expansion is computed Clm*cos + Slm*sin. In the serial case
       this gives me cache locality and thus more possibility to have both Clm
       and Slm for that l and m in cache.

        For the parallel case, each thread is going to compute Clm*cos and
       Slm*sin (or variants given by the Pines form of the Legendre expansion).
       This formatting guarantees locality of thread memory access.
    */

    // Be sure that buf isnt NULL
    if (buf != nullptr)
    {
        std::fstream file{path, std::ios::in};
        std::string buf_str;
        int i = 0;
        if (file.is_open())
        {
            // Begin reading file, stop at EOF or at degree
            while (std::getline(file, buf_str) && i < (bufsize - 2 * padding + 1))
            {
                // buf_str has now Clm,Slm
                auto start = 0U;
                auto end = buf_str.find(',');
                // Clm
                buf[i] = std::stod(buf_str.substr(start, end - start));

                start = end + 1;
                end = buf_str.size();
                // Slm
                buf[i + 1] = std::stod(buf_str.substr(start, end));
                i += 2;
            }
        }
        else
        { // File isnt good, throw FileError

            std::cout << "cof file not found" << std::endl;
            throw;
        }
    }
    else
    { // Pointer isnt good, throw BadPointer

        throw std::bad_alloc();
    }
}

size_t ceil2pow(uint16_t degree)
{
    if (degree <= 256)
    {
        return 256;
    }
    if (degree <= 512)
    {
        return 512;
    }
    if (degree <= 1024)
    {
        return 1024;
    }
    if (degree <= 2048)
    {
        return 2048;
    }
    return 0;
}

size_t ij2idx(int i, int j)
{
    return ((i) * (i + 1) / 2 + j);
}

// Load radius vectors from GMAT output file

void loadRadiusGMAT(double *&r, std::string path, size_t &N)
{
    /*
        GMAT output file doesnt give me the number of lines in the file
        itself, so im just eyeballing it and to avoid out of bounds
        access to r im just going to put a safeguard into the loop

    */

    // Open path
    std::fstream gmat_file{path};
    std::string buf_str;
    size_t i = 0;
    if (gmat_file.is_open())
    {
        // First line contains number of lines
        std::getline(gmat_file, buf_str);
        N = std::stoll(buf_str);
        r = new double[3 * N];

        while (std::getline(gmat_file, buf_str) && i < N)
        {
            // Separate string by space

            // Find first component
            auto start = 0U;
            auto end = buf_str.find(' ');
            r[i * 3] = std::stod(buf_str.substr(start, end - start));

            // Find second component

            start = end + 1;
            end = buf_str.find(' ', start);
            r[i * 3 + 1] = std::stod(buf_str.substr(start, end - start));

            // Now the last component
            start = end + 1;
            r[i * 3 + 2] = std::stod(buf_str.substr(start, buf_str.size()));
            i++;
        }
        if (i < N)
        {
            // The file has less lines than N is
            N = i;
        }
    }

    else
    {
        std::cout << "rad file not found" << std::endl;
        throw;
    }
}

void printMatrix(double *mat, uint16_t degree, std::string file)
{

    std::fstream f{file};
    for (int l = 0; l <= degree; l++)
    {
        for (int m = 0; m <= l; m++)
        {
            size_t idx = ij2idx(l, m);
            f << mat[idx] << ", " << std::endl;
        }
    }
}