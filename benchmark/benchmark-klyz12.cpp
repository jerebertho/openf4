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


int klyz12F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         KLYZ12                        #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    

    // Create variable name array                                         
    string * vars = new string[9];
    vars[0]='l';
    vars[1]='a';
    vars[2]='b';
    vars[3]='c';
    vars[4]='d';
    vars[5]='e';
    vars[6]='f';
    vars[7]='x';
    vars[8]='y';
    Monomial::setVariable(vars);

    // Init monomial tools
    Monomial::initMonomial(9);
    
    // Create polynomial array
    vector<Polynomial<eltType>> polKlyz12;
    
    // Fill the polynomial array
    polKlyz12.emplace_back("24*a^4*d^4+28*a^2*b^2*d^4+8*b^4*d^4+14*a^2*c^2*d^4+8*b^2*c^2*d^4+2*c^4*d^4+28*a^4*d^2*e^2+24*a^2*b^2*d^2*e^2+4*b^4*d^2*e^2+12*a^2*c^2*d^2*e^2+4*b^2*c^2*d^2*e^2+c^4*d^2*e^2+8*a^4*e^4+4*a^2*b^2*e^4+2*a^2*c^2*e^4+14*a^4*d^2*f^2+12*a^2*b^2*d^2*f^2+2*b^4*d^2*f^2+4*a^2*c^2*d^2*f^2+b^2*c^2*d^2*f^2+8*a^4*e^2*f^2+4*a^2*b^2*e^2*f^2+a^2*c^2*e^2*f^2+2*a^4*f^4+a^2*b^2*f^4+48*a^2*d^4*x+24*b^2*d^4*x+12*c^2*d^4*x+56*a^2*d^2*e^2*x+16*b^2*d^2*e^2*x+8*c^2*d^2*e^2*x+16*a^2*e^4*x+28*a^2*d^2*f^2*x+8*b^2*d^2*f^2*x+4*c^2*d^2*f^2*x+16*a^2*e^2*f^2*x+4*a^2*f^4*x+48*a^4*d^2*y+56*a^2*b^2*d^2*y+16*b^4*d^2*y+28*a^2*c^2*d^2*y+16*b^2*c^2*d^2*y+4*c^4*d^2*y+24*a^4*e^2*y+16*a^2*b^2*e^2*y+8*a^2*c^2*e^2*y+12*a^4*f^2*y+8*a^2*b^2*f^2*y+4*a^2*c^2*f^2*y+24*d^4*x^2+28*d^2*e^2*x^2+8*e^4*x^2+14*d^2*f^2*x^2+8*e^2*f^2*x^2+2*f^4*x^2+96*a^2*d^2*x*y+48*b^2*d^2*x*y+24*c^2*d^2*x*y+48*a^2*e^2*x*y+8*b^2*e^2*x*y+4*c^2*e^2*x*y+24*a^2*f^2*x*y+4*b^2*f^2*x*y+4*c^2*f^2*x*y+24*a^4*y^2+28*a^2*b^2*y^2+8*b^4*y^2+14*a^2*c^2*y^2+8*b^2*c^2*y^2+2*c^4*y^2+48*d^2*x^2*y+24*e^2*x^2*y+12*f^2*x^2*y+48*a^2*x*y^2+24*b^2*x*y^2+12*c^2*x*y^2+24*x^2*y^2");
    polKlyz12.emplace_back("96*l*a^3*d^4+56*l*a*b^2*d^4+28*l*a*c^2*d^4+112*l*a^3*d^2*e^2+48*l*a*b^2*d^2*e^2+24*l*a*c^2*d^2*e^2+32*l*a^3*e^4+8*l*a*b^2*e^4+4*l*a*c^2*e^4+56*l*a^3*d^2*f^2+24*l*a*b^2*d^2*f^2+8*l*a*c^2*d^2*f^2+32*l*a^3*e^2*f^2+8*l*a*b^2*e^2*f^2+2*l*a*c^2*e^2*f^2+8*l*a^3*f^4+2*l*a*b^2*f^4+96*l*a*d^4*x+112*l*a*d^2*e^2*x+32*l*a*e^4*x+56*l*a*d^2*f^2*x+32*l*a*e^2*f^2*x+8*l*a*f^4*x+192*l*a^3*d^2*y+112*l*a*b^2*d^2*y+56*l*a*c^2*d^2*y+96*l*a^3*e^2*y+32*l*a*b^2*e^2*y+16*l*a*c^2*e^2*y+48*l*a^3*f^2*y+16*l*a*b^2*f^2*y+8*l*a*c^2*f^2*y+192*l*a*d^2*x*y+96*l*a*e^2*x*y+48*l*a*f^2*x*y+96*l*a^3*y^2+56*l*a*b^2*y^2+28*l*a*c^2*y^2+96*l*a*x*y^2-1");
    polKlyz12.emplace_back("56*l*a^2*b*d^4+32*l*b^3*d^4+16*l*b*c^2*d^4+48*l*a^2*b*d^2*e^2+16*l*b^3*d^2*e^2+8*l*b*c^2*d^2*e^2+8*l*a^2*b*e^4+24*l*a^2*b*d^2*f^2+8*l*b^3*d^2*f^2+2*l*b*c^2*d^2*f^2+8*l*a^2*b*e^2*f^2+2*l*a^2*b*f^4+48*l*b*d^4*x+32*l*b*d^2*e^2*x+16*l*b*d^2*f^2*x+112*l*a^2*b*d^2*y+64*l*b^3*d^2*y+32*l*b*c^2*d^2*y+32*l*a^2*b*e^2*y+16*l*a^2*b*f^2*y+96*l*b*d^2*x*y+16*l*b*e^2*x*y+8*l*b*f^2*x*y+56*l*a^2*b*y^2+32*l*b^3*y^2+16*l*b*c^2*y^2+48*l*b*x*y^2-2");
    polKlyz12.emplace_back("28*l*a^2*c*d^4+16*l*b^2*c*d^4+8*l*c^3*d^4+24*l*a^2*c*d^2*e^2+8*l*b^2*c*d^2*e^2+4*l*c^3*d^2*e^2+4*l*a^2*c*e^4+8*l*a^2*c*d^2*f^2+2*l*b^2*c*d^2*f^2+2*l*a^2*c*e^2*f^2+24*l*c*d^4*x+16*l*c*d^2*e^2*x+8*l*c*d^2*f^2*x+56*l*a^2*c*d^2*y+32*l*b^2*c*d^2*y+16*l*c^3*d^2*y+16*l*a^2*c*e^2*y+8*l*a^2*c*f^2*y+48*l*c*d^2*x*y+8*l*c*e^2*x*y+8*l*c*f^2*x*y+28*l*a^2*c*y^2+16*l*b^2*c*y^2+8*l*c^3*y^2+24*l*c*x*y^2-3");
    polKlyz12.emplace_back("96*l*a^4*d^3+112*l*a^2*b^2*d^3+32*l*b^4*d^3+56*l*a^2*c^2*d^3+32*l*b^2*c^2*d^3+8*l*c^4*d^3+56*l*a^4*d*e^2+48*l*a^2*b^2*d*e^2+8*l*b^4*d*e^2+24*l*a^2*c^2*d*e^2+8*l*b^2*c^2*d*e^2+2*l*c^4*d*e^2+28*l*a^4*d*f^2+24*l*a^2*b^2*d*f^2+4*l*b^4*d*f^2+8*l*a^2*c^2*d*f^2+2*l*b^2*c^2*d*f^2+192*l*a^2*d^3*x+96*l*b^2*d^3*x+48*l*c^2*d^3*x+112*l*a^2*d*e^2*x+32*l*b^2*d*e^2*x+16*l*c^2*d*e^2*x+56*l*a^2*d*f^2*x+16*l*b^2*d*f^2*x+8*l*c^2*d*f^2*x+96*l*a^4*d*y+112*l*a^2*b^2*d*y+32*l*b^4*d*y+56*l*a^2*c^2*d*y+32*l*b^2*c^2*d*y+8*l*c^4*d*y+96*l*d^3*x^2+56*l*d*e^2*x^2+28*l*d*f^2*x^2+192*l*a^2*d*x*y+96*l*b^2*d*x*y+48*l*c^2*d*x*y+96*l*d*x^2*y-4");
    polKlyz12.emplace_back("56*l*a^4*d^2*e+48*l*a^2*b^2*d^2*e+8*l*b^4*d^2*e+24*l*a^2*c^2*d^2*e+8*l*b^2*c^2*d^2*e+2*l*c^4*d^2*e+32*l*a^4*e^3+16*l*a^2*b^2*e^3+8*l*a^2*c^2*e^3+16*l*a^4*e*f^2+8*l*a^2*b^2*e*f^2+2*l*a^2*c^2*e*f^2+112*l*a^2*d^2*e*x+32*l*b^2*d^2*e*x+16*l*c^2*d^2*e*x+64*l*a^2*e^3*x+32*l*a^2*e*f^2*x+48*l*a^4*e*y+32*l*a^2*b^2*e*y+16*l*a^2*c^2*e*y+56*l*d^2*e*x^2+32*l*e^3*x^2+16*l*e*f^2*x^2+96*l*a^2*e*x*y+16*l*b^2*e*x*y+8*l*c^2*e*x*y+48*l*e*x^2*y-5");
    polKlyz12.emplace_back("28*l*a^4*d^2*f+24*l*a^2*b^2*d^2*f+4*l*b^4*d^2*f+8*l*a^2*c^2*d^2*f+2*l*b^2*c^2*d^2*f+16*l*a^4*e^2*f+8*l*a^2*b^2*e^2*f+2*l*a^2*c^2*e^2*f+8*l*a^4*f^3+4*l*a^2*b^2*f^3+56*l*a^2*d^2*f*x+16*l*b^2*d^2*f*x+8*l*c^2*d^2*f*x+32*l*a^2*e^2*f*x+16*l*a^2*f^3*x+24*l*a^4*f*y+16*l*a^2*b^2*f*y+8*l*a^2*c^2*f*y+28*l*d^2*f*x^2+16*l*e^2*f*x^2+8*l*f^3*x^2+48*l*a^2*f*x*y+8*l*b^2*f*x*y+8*l*c^2*f*x*y+24*l*f*x^2*y-6");
    polKlyz12.emplace_back("48*l*a^2*d^4+24*l*b^2*d^4+12*l*c^2*d^4+56*l*a^2*d^2*e^2+16*l*b^2*d^2*e^2+8*l*c^2*d^2*e^2+16*l*a^2*e^4+28*l*a^2*d^2*f^2+8*l*b^2*d^2*f^2+4*l*c^2*d^2*f^2+16*l*a^2*e^2*f^2+4*l*a^2*f^4+48*l*d^4*x+56*l*d^2*e^2*x+16*l*e^4*x+28*l*d^2*f^2*x+16*l*e^2*f^2*x+4*l*f^4*x+96*l*a^2*d^2*y+48*l*b^2*d^2*y+24*l*c^2*d^2*y+48*l*a^2*e^2*y+8*l*b^2*e^2*y+4*l*c^2*e^2*y+24*l*a^2*f^2*y+4*l*b^2*f^2*y+4*l*c^2*f^2*y+96*l*d^2*x*y+48*l*e^2*x*y+24*l*f^2*x*y+48*l*a^2*y^2+24*l*b^2*y^2+12*l*c^2*y^2+48*l*x*y^2-7");
    polKlyz12.emplace_back("48*l*a^4*d^2+56*l*a^2*b^2*d^2+16*l*b^4*d^2+28*l*a^2*c^2*d^2+16*l*b^2*c^2*d^2+4*l*c^4*d^2+24*l*a^4*e^2+16*l*a^2*b^2*e^2+8*l*a^2*c^2*e^2+12*l*a^4*f^2+8*l*a^2*b^2*f^2+4*l*a^2*c^2*f^2+96*l*a^2*d^2*x+48*l*b^2*d^2*x+24*l*c^2*d^2*x+48*l*a^2*e^2*x+8*l*b^2*e^2*x+4*l*c^2*e^2*x+24*l*a^2*f^2*x+4*l*b^2*f^2*x+4*l*c^2*f^2*x+48*l*a^4*y+56*l*a^2*b^2*y+16*l*b^4*y+28*l*a^2*c^2*y+16*l*b^2*c^2*y+4*l*c^4*y+48*l*d^2*x^2+24*l*e^2*x^2+12*l*f^2*x^2+96*l*a^2*x*y+48*l*b^2*x*y+24*l*c^2*x*y+48*l*x^2*y-8");

    // Create klyz12 ideal;
    Ideal<eltType> klyz12(polKlyz12, 9, 100000);
    
    // Compute a reduced groebner basis;
    nbGen=klyz12.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        klyz12.printReducedGroebnerBasis("benchmark-klyz12.res", modulo);
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
    nbGen=klyz12F4(magma);

    cout << "klyz12: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


