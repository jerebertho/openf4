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


int kat9F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         KAT9                          #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(9);
    
    // Create variable name array                                                                              
    string * vars = new string[9];
    for(int i = 0; i < 9; i++)
      {
        vars[i]='x'+to_string(i+1);
      }
    Monomial::setVariable(vars);

    // Create polynomial array
    vector<Polynomial<eltType>> polKat9;
    
    // Fill the polynomial array
    polKat9.emplace_back("x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9-1");
    polKat9.emplace_back("x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2-x1");
    polKat9.emplace_back("2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9-x2");
    polKat9.emplace_back("x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9-x3");
    polKat9.emplace_back("2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9-x4");
    polKat9.emplace_back("x3^2+2*x2*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9-x5");
    polKat9.emplace_back("2*x3*x4+2*x2*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9-x6");
    polKat9.emplace_back("x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x2*x8+2*x3*x9-x7");
    polKat9.emplace_back("2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x2*x9-x8");

    // Create kat9 ideal;
    Ideal<eltType> kat9(polKat9, 9 ,1000000);
    
    // Compute a reduced groebner basis;
    nbGen=kat9.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        kat9.printReducedGroebnerBasis("benchmark-katsura9", modulo);
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
    nbGen=kat9F4(magma);

    cout << "kat9: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


