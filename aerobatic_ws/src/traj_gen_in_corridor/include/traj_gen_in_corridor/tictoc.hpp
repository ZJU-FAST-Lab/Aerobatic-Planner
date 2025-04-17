#ifndef TICTOC_HPP
#define TICTOC_HPP

#include <ctime>
#include <cstdlib>
#include <chrono>

// Timer for computation efficiency profile
class TicToc
{
public:
    inline void tic()
    {
        ini = std::chrono::high_resolution_clock::now();
    }

    inline double toc()
    {
        fin = std::chrono::high_resolution_clock::now();
        return 1000.0 * std::chrono::duration_cast<std::chrono::duration<double>>(fin - ini).count();
    }

    TicToc()
    {
        tic();
    }

private:
    std::chrono::high_resolution_clock::time_point ini, fin;
};

#endif