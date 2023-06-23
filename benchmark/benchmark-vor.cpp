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


int vorF4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         VOR                           #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    

    // Create variable name array                                         
    string * vars = new string[8];
    for(int i = 0; i < 8; i++)
      {
	vars[i]='z'+to_string(i+1);
      }
    Monomial::setVariable(vars);

    // Init monomial tools
    Monomial::initMonomial(8);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polVor;
    
    // Fill the polynomial array
    polVor.emplace_back("64*z2^4*z3^2*z5^2-256*z2^3*z3*z4*z5^3+256*z2^2*z3^2*z5^4+256*z2^2*z4^2*z5^4+64*z2^4*z3^2*z5-384*z2^3*z3*z4*z5^2+512*z2^2*z3^2*z5^3+512*z2^2*z4^2*z5^3+64*z2^4*z3*z5*z6-128*z2^3*z4*z5^2*z6+256*z2^2*z3*z5^3*z6-128*z2^3*z3*z5^2*z7+256*z2^2*z4*z5^3*z7+16*z2^4*z3^2-192*z2^3*z3*z4*z5+384*z2^2*z3^2*z5^2+384*z2^2*z4^2*z5^2-256*z2*z3*z4*z5^3+256*z2^2*z5^4+32*z2^4*z3*z6-128*z2^3*z4*z5*z6+384*z2^2*z3*z5^2*z6+16*z2^4*z6^2+64*z2^2*z5^2*z6^2-128*z2^3*z3*z5*z7+384*z2^2*z4*z5^2*z7-64*z2^3*z5*z6*z7+64*z2^2*z5^2*z7^2-32*z2^3*z3*z4+128*z2^2*z3^2*z5+128*z2^2*z4^2*z5-384*z2*z3*z4*z5^2+512*z2^2*z5^3-32*z2^3*z4*z6+192*z2^2*z3*z5*z6-128*z2*z4*z5^2*z6+64*z2^2*z5*z6^2-32*z2^3*z3*z7+192*z2^2*z4*z5*z7-128*z2*z3*z5^2*z7-32*z2^3*z6*z7+64*z2^2*z5*z7^2+16*z2^2*z3^2+16*z2^2*z4^2-192*z2*z3*z4*z5+256*z2^2*z5^2+64*z4^2*z5^2+32*z2^2*z3*z6-128*z2*z4*z5*z6+16*z2^2*z6^2+32*z2^2*z4*z7-128*z2*z3*z5*z7-64*z2*z5*z6*z7+16*z2^2*z7^2-32*z2*z3*z4+64*z4^2*z5-32*z2*z4*z6-32*z2*z3*z7+64*z4*z5*z7-32*z2*z6*z7+16*z4^2+32*z4*z7+16*z7^2+1");
    polVor.emplace_back("16*z1*z2^3*z3^2*z5^2-48*z1*z2^2*z3*z4*z5^3+32*z1*z2*z3^2*z5^4+32*z1*z2*z4^2*z5^4+16*z1*z2^3*z3^2*z5-72*z1*z2^2*z3*z4*z5^2+64*z1*z2*z3^2*z5^3+64*z1*z2*z4^2*z5^3+16*z1*z2^3*z3*z5*z6-24*z1*z2^2*z4*z5^2*z6+32*z1*z2*z3*z5^3*z6-24*z1*z2^2*z3*z5^2*z7+32*z1*z2*z4*z5^3*z7+4*z1*z2^3*z3^2-36*z1*z2^2*z3*z4*z5+48*z1*z2*z3^2*z5^2+48*z1*z2*z4^2*z5^2-16*z1*z3*z4*z5^3+32*z1*z2*z5^4+8*z1*z2^3*z3*z6-24*z1*z2^2*z4*z5*z6+48*z1*z2*z3*z5^2*z6+4*z1*z2^3*z6^2+8*z1*z2*z5^2*z6^2-24*z1*z2^2*z3*z5*z7+48*z1*z2*z4*z5^2*z7-12*z1*z2^2*z5*z6*z7+8*z1*z2*z5^2*z7^2-6*z1*z2^2*z3*z4+16*z1*z2*z3^2*z5+16*z1*z2*z4^2*z5-24*z1*z3*z4*z5^2+64*z1*z2*z5^3-6*z1*z2^2*z4*z6+24*z1*z2*z3*z5*z6-8*z1*z4*z5^2*z6+8*z1*z2*z5*z6^2-6*z1*z2^2*z3*z7+24*z1*z2*z4*z5*z7-8*z1*z3*z5^2*z7-6*z1*z2^2*z6*z7+8*z1*z2*z5*z7^2+2*z1*z2*z3^2+2*z1*z2*z4^2-12*z1*z3*z4*z5+32*z1*z2*z5^2+4*z1*z2*z3*z6-8*z1*z4*z5*z6+2*z1*z2*z6^2+4*z1*z2*z4*z7-8*z1*z3*z5*z7-4*z1*z5*z6*z7+2*z1*z2*z7^2-2*z1*z3*z4-2*z1*z4*z6-2*z1*z3*z7-2*z1*z6*z7-z2+1");
    polVor.emplace_back("8*z1*z2^4*z3*z5^2-16*z1*z2^3*z4*z5^3+32*z1*z2^2*z3*z5^4+8*z1*z2^4*z3*z5-24*z1*z2^3*z4*z5^2+64*z1*z2^2*z3*z5^3+4*z1*z2^4*z5*z6+16*z1*z2^2*z5^3*z6-8*z1*z2^3*z5^2*z7+2*z1*z2^4*z3-12*z1*z2^3*z4*z5+48*z1*z2^2*z3*z5^2-16*z1*z2*z4*z5^3+2*z1*z2^4*z6+24*z1*z2^2*z5^2*z6-8*z1*z2^3*z5*z7-2*z1*z2^3*z4+16*z1*z2^2*z3*z5-24*z1*z2*z4*z5^2+12*z1*z2^2*z5*z6-2*z1*z2^3*z7-8*z1*z2*z5^2*z7+2*z1*z2^2*z3-12*z1*z2*z4*z5+2*z1*z2^2*z6-8*z1*z2*z5*z7-2*z1*z2*z4-2*z1*z2*z7-z3+2");
    polVor.emplace_back("-16*z1*z2^3*z3*z5^3+32*z1*z2^2*z4*z5^4-24*z1*z2^3*z3*z5^2+64*z1*z2^2*z4*z5^3-8*z1*z2^3*z5^2*z6+16*z1*z2^2*z5^3*z7-12*z1*z2^3*z3*z5+48*z1*z2^2*z4*z5^2-16*z1*z2*z3*z5^3-8*z1*z2^3*z5*z6+24*z1*z2^2*z5^2*z7-2*z1*z2^3*z3+16*z1*z2^2*z4*z5-24*z1*z2*z3*z5^2-2*z1*z2^3*z6-8*z1*z2*z5^2*z6+12*z1*z2^2*z5*z7+2*z1*z2^2*z4-12*z1*z2*z3*z5+8*z1*z4*z5^2-8*z1*z2*z5*z6+2*z1*z2^2*z7-2*z1*z2*z3+8*z1*z4*z5-2*z1*z2*z6+4*z1*z5*z7+2*z1*z4+2*z1*z7-z4+6");
    polVor.emplace_back("8*z1*z2^4*z3^2*z5-48*z1*z2^3*z3*z4*z5^2+64*z1*z2^2*z3^2*z5^3+64*z1*z2^2*z4^2*z5^3+4*z1*z2^4*z3^2-48*z1*z2^3*z3*z4*z5+96*z1*z2^2*z3^2*z5^2+96*z1*z2^2*z4^2*z5^2+4*z1*z2^4*z3*z6-16*z1*z2^3*z4*z5*z6+48*z1*z2^2*z3*z5^2*z6-16*z1*z2^3*z3*z5*z7+48*z1*z2^2*z4*z5^2*z7-12*z1*z2^3*z3*z4+48*z1*z2^2*z3^2*z5+48*z1*z2^2*z4^2*z5-48*z1*z2*z3*z4*z5^2+64*z1*z2^2*z5^3-8*z1*z2^3*z4*z6+48*z1*z2^2*z3*z5*z6+8*z1*z2^2*z5*z6^2-8*z1*z2^3*z3*z7+48*z1*z2^2*z4*z5*z7-4*z1*z2^3*z6*z7+8*z1*z2^2*z5*z7^2+8*z1*z2^2*z3^2+8*z1*z2^2*z4^2-48*z1*z2*z3*z4*z5+96*z1*z2^2*z5^2+12*z1*z2^2*z3*z6-16*z1*z2*z4*z5*z6+4*z1*z2^2*z6^2+12*z1*z2^2*z4*z7-16*z1*z2*z3*z5*z7+4*z1*z2^2*z7^2-12*z1*z2*z3*z4+32*z1*z2^2*z5+8*z1*z4^2*z5-8*z1*z2*z4*z6-8*z1*z2*z3*z7-4*z1*z2*z6*z7+4*z1*z4^2+4*z1*z4*z7-z5+24");
    polVor.emplace_back("4*z1*z2^4*z3*z5-8*z1*z2^3*z4*z5^2+16*z1*z2^2*z3*z5^3+2*z1*z2^4*z3-8*z1*z2^3*z4*z5+24*z1*z2^2*z3*z5^2+2*z1*z2^4*z6+8*z1*z2^2*z5^2*z6-4*z1*z2^3*z5*z7-2*z1*z2^3*z4+12*z1*z2^2*z3*z5-8*z1*z2*z4*z5^2+8*z1*z2^2*z5*z6-2*z1*z2^3*z7+2*z1*z2^2*z3-8*z1*z2*z4*z5+2*z1*z2^2*z6-4*z1*z2*z5*z7-2*z1*z2*z4-2*z1*z2*z7-z6+120");
    polVor.emplace_back("-8*z1*z2^3*z3*z5^2+16*z1*z2^2*z4*z5^3-8*z1*z2^3*z3*z5+24*z1*z2^2*z4*z5^2-4*z1*z2^3*z5*z6+8*z1*z2^2*z5^2*z7-2*z1*z2^3*z3+12*z1*z2^2*z4*z5-8*z1*z2*z3*z5^2-2*z1*z2^3*z6+8*z1*z2^2*z5*z7+2*z1*z2^2*z4-8*z1*z2*z3*z5-4*z1*z2*z5*z6+2*z1*z2^2*z7-2*z1*z2*z3+4*z1*z4*z5-2*z1*z2*z6+2*z1*z4+2*z1*z7-z7+720");
    polVor.emplace_back("z1+z2+z3+z4+z5+z6+z7+z8");

    // Create vor ideal;
    Ideal<eltType> vor(polVor, 8, 100000);
    
    // Compute a reduced groebner basis;
    nbGen=vor.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        vor.printReducedGroebnerBasis("benchmark-vor.res", modulo);
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
    nbGen=vorF4(magma);

    cout << "vor: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


