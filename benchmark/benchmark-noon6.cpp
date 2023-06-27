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


int noon6F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         NOON6                         #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    

    // Create variable name array                                         
    string * vars = new string[6];
    for(int i = 0; i < 6; i++)
      {
        vars[i]='x'+to_string(i+1);
      }
    Monomial::setVariable(vars);

    // Init monomial tools
    Monomial::initMonomial(8);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polNoon6;
    
    // Fill the polynomial array
    polNoon6.emplace_back("10*x1^2*x6+10*x2^2*x6+10*x3^2*x6+10*x4^2*x6+10*x5^2*x6-11*x6+10");
    polNoon6.emplace_back("10*x1^2*x5+10*x2^2*x5+10*x3^2*x5+10*x4^2*x5+10*x5*x6^2-11*x5+10");
    polNoon6.emplace_back("10*x1^2*x4+10*x2^2*x4+10*x3^2*x4+10*x4*x5^2+10*x4*x6^2-11*x4+10");
    polNoon6.emplace_back("10*x1^2*x3+10*x2^2*x3+10*x3*x4^2+10*x3*x5^2+10*x3*x6^2-11*x3+10");
    polNoon6.emplace_back("10*x1*x2^2+10*x1*x3^2+10*x1*x4^2+10*x1*x5^2+10*x1*x6^2-11*x1+10");
    polNoon6.emplace_back("10*x1^2*x2+10*x2*x3^2+10*x2*x4^2+10*x2*x5^2+10*x2*x6^2-11*x2+10");

    // Create noon6 ideal;
    Ideal<eltType> noon6(polNoon6, 8, 100000);
    
    // Compute a reduced groebner basis;
    nbGen=noon6.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        noon6.printReducedGroebnerBasis("benchmark-noon6.res", modulo);
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
    nbGen=noon6F4(magma);

    cout << "noon6: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


