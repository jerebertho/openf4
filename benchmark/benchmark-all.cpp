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
 *  \brief Benchmark with integer 16 bits coefficients.
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
typedef ElementPrime<int16_t> eltType;
int16_t modulo=3;


int allF4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         ALL                           #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(12);
    
    // Create variable name array                                         
    string * vars = new string[12];
    for(int i = 0; i < 12; i++)
      {
	vars[i]='t'+to_string(i+1);
      }
    Monomial::setVariable(vars);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polAll;
    
    // Fill the polynomial array
    polAll.emplace_back("t10^3*t11 - t10*t11 - 1");
    polAll.emplace_back("-t2*t4^3*t9^2*t10^2*t12^6 + t2*t4*t9^4*t10^4*t12^6 - t3^2*t9^4*t10^4*t12^6 + t3*t4^3*t8*t9*t10^2*t12^6 + t3*t4*t8*t9^3*t10^4*t12^6 - t4^4*t8^2*t10^2*t12^6 + 1");
    polAll.emplace_back("-t1*t2*t9^4*t10^2 - t1*t3*t8*t9^3*t10^2 + t1*t4^3*t8^2 + t2^2*t8*t9^3*t10^2 + t2*t3^3*t9^2 - t2*t4^3*t6*t9 - t2*t4*t6*t9^3*t10^2 - t2*t4*t8^3*t9*t10^2 - t3^4*t8*t9 + t3^3*t4*t8^2 + t3^2*t6*t9^3*t10^2 + t3^2*t8^3*t9*t10^2 - t3*t4^3*t6*t8 - t3*t4*t8^4*t10^2 - t4^4*t9 + t4^3*t9^2 + t4^2*t9^3*t10^2 - t4*t9^4*t10^2");
    polAll.emplace_back("t1*t2*t6*t9^3*t10^2 + t1*t2*t8^3*t9*t10^2 - t1*t3^3*t8^2 + t1*t3*t8^4*t10^2 + t1*t4^3*t9 + t1*t4*t9^3*t10^2 + t1*t9^4*t10^2 - t2^4*t9^2 + t2^3*t3*t8*t9 - t2^3*t4*t8^2 - t2^2*t8^4*t10^2 + t2*t3^3*t6*t9 + t2*t3*t9^3*t10^2 - t2*t4^3*t6^2 + t2*t4^3*t8 + t2*t4*t6*t8^3*t10^2 + t2*t8*t9^3*t10^2 + t3^4*t6*t8 + t3^3*t4*t9 - t3^3*t9^2 - t3^2*t6*t8^3*t10^2 + t4^4*t6 + t4^3*t6*t9 - t4^2*t8^3*t10^2 + t4*t6*t9^3*t10^2 + t4*t8^3*t9*t10^2");
    polAll.emplace_back("t1^3*t2*t9^2 - t1^3*t3*t8*t9 + t1^3*t4*t8^2 + t1^2*t9^3*t10^2 + t1*t2^3*t8^2 - t1*t2*t6*t8^3*t10^2 - t1*t3^3*t9 - t1*t4^3*t6 - t1*t4*t8^3*t10^2 - t1*t6*t9^3*t10^2 - t1*t8^3*t9*t10^2 - t2^4*t6*t9 - t2^3*t3*t6*t8 - t2^3*t4*t9 + t2^3*t9^2 + t2*t3^3*t6^2 - t2*t3^3*t8 - t2*t3*t8^3*t10^2 - t2*t4*t6^3*t9*t10^2 - t2*t8^4*t10^2 - t3^3*t4*t6 - t3^3*t6*t9 + t3^2*t6^3*t9*t10^2 + t3*t4^3 - t3*t4*t6^3*t8*t10^2 - t3*t9^3*t10^2 + t4^3*t6^2 - t4^3*t8 - t4*t6*t8^3*t10^2 + t8*t9^3*t10^2");
    polAll.emplace_back("-t1^4*t8^2 + t1^3*t2*t6*t9 + t1^3*t3*t6*t8 + t1^3*t4*t9 - t1^3*t9^2 - t1^2*t8^3*t10^2 + t1*t2^3*t9 + t1*t2*t6^3*t9*t10^2 + t1*t3^3*t6 + t1*t3*t6^3*t8*t10^2 + t1*t6*t8^3*t10^2 - t2^4*t6^2 + t2^4*t8 + t2^3*t4*t6 + t2^3*t6*t9 - t2^2*t6^3*t8*t10^2 + t2*t4*t6^4*t10^2 - t3^4 - t3^3*t6^2 + t3^3*t8 - t3^2*t6^4*t10^2 + t3*t8^3*t10^2 - t4^2*t6^3*t10^2 + t4*t6^3*t9*t10^2 - t8^4*t10^2");
    polAll.emplace_back("-t1^4*t9 + t1^3*t2*t6^2 - t1^3*t2*t8 - t1^3*t4*t6 - t1^3*t6*t9 - t1*t2^3*t6 - t1*t2*t6^4*t10^2 - t1*t4*t6^3*t10^2 - t1*t6^3*t9*t10^2 + t2^3*t3 + t2^3*t6^2 - t2^3*t8 - t2*t3*t6^3*t10^2 - t2*t4*t9*t10^2 - t2*t6^3*t8*t10^2 + t2*t9^2 + t3^2*t9*t10^2 - t3*t4*t8*t10^2 - t3*t8*t9 - t4*t6^4*t10^2 + t4*t8^2");
    polAll.emplace_back("t1^4*t6 - t1^3*t3 - t1^3*t6^2 + t1^3*t8 - t1^2*t6^3*t10^2 + t1*t2*t9*t10^2 + t1*t3*t8*t10^2 + t1*t6^4*t10^2 - t1*t8^2 - t2^2*t8*t10^2 + t2*t4*t6*t10^2 + t2*t6*t9 - t3^2*t6*t10^2 + t3*t6^3*t10^2 + t3*t6*t8 - t4^2*t10^2 + t4*t9*t10^2 + t4*t9 - t6^3*t8*t10^2 - t9^2");
    polAll.emplace_back("-t1*t2*t6*t10^2 - t1*t4*t10^2 - t1*t9*t10^2 - t1*t9 - t2*t3*t10^2 + t2*t6^2 - t2*t8*t10^2 - t2*t8 - t4*t6*t10^2 - t4*t6 - t6*t9");
    polAll.emplace_back("-t1^2*t10^2 + t1*t6*t10^2 + t1*t6 + t3*t10^2 - t3 - t6^2 - t8*t10^2 + t8");

    // Create All ideal;
    Ideal<eltType> all(polAll, 12, 10000000);
    
    // Compute a reduced groebner basis;
    nbGen=all.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        all.printReducedGroebnerBasis("benchmark-all.res", modulo);
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
    nbGen=allF4(magma);

    cout << "all: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


