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



int eco11F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         ECO11                         #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(11);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polEco11;
    
    // Fill the polynomial array
    polEco11.emplace_back("x0*x1*x10+x1*x2*x10+x2*x3*x10+x3*x4*x10+x4*x5*x10+x5*x6*x10+x6*x7*x10+x7*x8*x10+x8*x9*x10+x0*x10-1");
    polEco11.emplace_back("x0*x2*x10+x1*x3*x10+x2*x4*x10+x3*x5*x10+x4*x6*x10+x5*x7*x10+x6*x8*x10+x7*x9*x10+x1*x10-2");
    polEco11.emplace_back("x0*x3*x10+x1*x4*x10+x2*x5*x10+x3*x6*x10+x4*x7*x10+x5*x8*x10+x6*x9*x10+x2*x10-3");
    polEco11.emplace_back("x0*x4*x10+x1*x5*x10+x2*x6*x10+x3*x7*x10+x4*x8*x10+x5*x9*x10+x3*x10-4");
    polEco11.emplace_back("x0*x5*x10+x1*x6*x10+x2*x7*x10+x3*x8*x10+x4*x9*x10+x4*x10-5");
    polEco11.emplace_back("x0*x6*x10+x1*x7*x10+x2*x8*x10+x3*x9*x10+x5*x10-6");
    polEco11.emplace_back("x0*x7*x10+x1*x8*x10+x2*x9*x10+x6*x10-7");
    polEco11.emplace_back("x0*x8*x10+x1*x9*x10+x7*x10-8");
    polEco11.emplace_back("x0*x9*x10+x8*x10-9");
    polEco11.emplace_back("x9*x10-10");
    polEco11.emplace_back("x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+1");

    // Create eco11 ideal;
    Ideal<eltType> eco11(polEco11, 11, 11000000);
    
    // Compute a reduced groebner basis;
    nbGen=eco11.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        eco11.printReducedGroebnerBasis("benchmark-eco11.res", modulo);
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
    nbGen=eco11F4(magma);

    cout << "eco11: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


