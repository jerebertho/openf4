/* 
 * Copyright (C) 2015 Antoine Joux, Vanessa Vitse and Titouan Coladon
 * 
 * This file is part of openf4.
 * 
 * openf4 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * openf4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with openf4.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 *  \file benchmark-long.cpp
 *  \example benchmark-long.cpp
 *  \brief Benchmark with integer 64 bits coefficients.
 *  \ingroup benchmark
 *  \author Vanessa VITSE, Antoine JOUX, Titouan COLADON
 */

#include <iostream>
#include <openf4.h>

using namespace F4;
using namespace std;

// Global variable
int F4::VERBOSE=1;
#ifdef USE_OPENMP
int F4::NB_THREAD=min(1/*8*/, omp_get_num_procs());
#else
int F4::NB_THREAD=1;
#endif

// Init element-prime tools
typedef ElementPrime<int64_t> eltType;
int64_t modulo=1073741827LL;



int eco12F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         ECO12                         #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(12);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polEco12;
    
    // Fill the polynomial array
    polEco12.emplace_back("x0*x1*x11+x1*x2*x11+x2*x3*x11+x3*x4*x11+x4*x5*x11+x5*x6*x11+x6*x7*x11+x7*x8*x11+x8*x9*x11+x9*x10*x11+x0*x11-1");
    polEco12.emplace_back("x0*x2*x11+x1*x3*x11+x2*x4*x11+x3*x5*x11+x4*x6*x11+x5*x7*x11+x6*x8*x11+x7*x9*x11+x8*x10*x11+x1*x11-2");
    polEco12.emplace_back("x0*x3*x11+x1*x4*x11+x2*x5*x11+x3*x6*x11+x4*x7*x11+x5*x8*x11+x6*x9*x11+x7*x10*x11+x2*x11-3");
    polEco12.emplace_back("x0*x4*x11+x1*x5*x11+x2*x6*x11+x3*x7*x11+x4*x8*x11+x5*x9*x11+x6*x10*x11+x3*x11-4");
    polEco12.emplace_back("x0*x5*x11+x1*x6*x11+x2*x7*x11+x3*x8*x11+x4*x9*x11+x5*x10*x11+x4*x11-5");
    polEco12.emplace_back("x0*x6*x11+x1*x7*x11+x2*x8*x11+x3*x9*x11+x4*x10*x11+x5*x11-6");
    polEco12.emplace_back("x0*x7*x11+x1*x8*x11+x2*x9*x11+x3*x10*x11+x6*x11-7");
    polEco12.emplace_back("x0*x8*x11+x1*x9*x11+x2*x10*x11+x7*x11-8");
    polEco12.emplace_back("x0*x9*x11+x1*x10*x11+x8*x11-9");
    polEco12.emplace_back("x0*x10*x11+x9*x11-10");
    polEco12.emplace_back("x10*x11-11");
    polEco12.emplace_back("x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+1");

    // Create eco12 ideal;
    Ideal<eltType> eco12(polEco12, 12, 12000000);
    
    // Compute a reduced groebner basis;
    nbGen=eco12.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        eco12.printReducedGroebnerBasis("benchmark-eco12.res", modulo);
    }
    return nbGen;
}

int main (int argc, char **argv)
{
    // Time
    chrono::steady_clock::time_point start;
    typedef chrono::duration<int,milli> millisecs_t;
    
    // Number of thread(s)
    cout << NB_THREAD << " thread(s) used " << endl << endl;
    
    // Magma output
    bool magma = false;
    
    // Number of generator
    int nbGen;
    
    cout << "Benchmark for ideal with integer long type coefficient." << endl;
    
    start=chrono::steady_clock::now();
    nbGen=eco12F4(magma);

    cout << "eco12: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


