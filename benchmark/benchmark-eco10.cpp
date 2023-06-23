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



int eco10F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         ECO10                         #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(10);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polEco10;
    
    // Fill the polynomial array
    polEco10.emplace_back("x0*x1*x9+x1*x2*x9+x2*x3*x9+x3*x4*x9+x4*x5*x9+x5*x6*x9+x6*x7*x9+x7*x8*x9+x0*x9-1");
    polEco10.emplace_back("x0*x2*x9+x1*x3*x9+x2*x4*x9+x3*x5*x9+x4*x6*x9+x5*x7*x9+x6*x8*x9+x1*x9-2");
    polEco10.emplace_back("x0*x3*x9+x1*x4*x9+x2*x5*x9+x3*x6*x9+x4*x7*x9+x5*x8*x9+x2*x9-3");
    polEco10.emplace_back("x0*x4*x9+x1*x5*x9+x2*x6*x9+x3*x7*x9+x4*x8*x9+x3*x9-4");
    polEco10.emplace_back("x0*x5*x9+x1*x6*x9+x2*x7*x9+x3*x8*x9+x4*x9-5");
    polEco10.emplace_back("x0*x6*x9+x1*x7*x9+x2*x8*x9+x5*x9-6");
    polEco10.emplace_back("x0*x7*x9+x1*x8*x9+x6*x9-7");
    polEco10.emplace_back("x0*x8*x9+x7*x9-8");
    polEco10.emplace_back("x8*x9-9");
    polEco10.emplace_back("x0+x1+x2+x3+x4+x5+x6+x7+x8+1");

    // Create eco10 ideal;
    Ideal<eltType> eco10(polEco10, 10, 10000000);
    
    // Compute a reduced groebner basis;
    nbGen=eco10.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        eco10.printReducedGroebnerBasis("benchmark-eco10.res", modulo);
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
    nbGen=eco10F4(magma);

    cout << "eco10: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


