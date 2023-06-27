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


int cim3RF4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         CIM3R                         #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(5);
    
    // Create variable name array                                                                              
    string * vars = new string[5];
    vars[0]='x';
    vars[1]='y';
    vars[2]='z';
    vars[3]='w';
    vars[4]='t';
    Monomial::setVariable(vars);

    // Create polynomial array
    vector<Polynomial<eltType>> polCim3R;
    
    // Fill the polynomial array
    polCim3R.emplace_back("-10402028*x^4 + 352976875*x^3 - 471760098*x^2 - 503625715*x - 87412*y^4 - 188275063*y^3 - 536075092*y^2 + 177776532*y + 10402028*z^4 - 326739522*z^3 - 242162594*z^2 + 447033182*z + 87412*w^4 - 226530569*w^3 - 467660625*w^2 + 66858318*w + 136762709");
    polCim3R.emplace_back("87412*x^4 - 327795641*x^3 - 338910092*x^2 - 169426269*x + 197223747*y^8 - 535572841*y^7 - 215411194*y^6 + 440229141*y^5 + 438678440*y^4 + 506929309*y^3 + 97902739*y^2 - 468625873*y - 87412*z^4 + 29814831*z^3 - 223542093*z^2 - 413081344*z - 197223747*w^8 + 74203091*w^7 + 87409870*w^6 + 126952622*w^5 + 242460333*w^4 - 398382159*w^3 + 231713254*w^2 - 41203551*w + 27014171");
    polCim3R.emplace_back("-197223747*x^8 + 185160619*x^7 - 493773767*x^6 + 170789484*x^5 + 311453025*x^4 + 446564908*x^3 + 346516772*x^2 + 244751246*x + 275991420*y^12 + 526765638*y^11 - 477444453*y^10 + 148838230*y^9 - 105083472*y^8 + 251903780*y^7 - 456007032*y^6 - 330333080*y^5 + 215329001*y^4 + 319999046*y^3 - 491268447*y^2 + 114424000*y - 192202996");
    polCim3R.emplace_back("-197223747*z^8 + 79238883*z^7 + 203714910*z^6 - 346976752*z^5 - 188187845*z^4 + 526426098*z^3 + 21501468*z^2 + 400413782*z + 275991420*w^12 + 186535609*w^11 + 280333695*w^10 + 241314867*w^9 - 329065829*w^8 + 32937494*w^7 + 407769360*w^6 + 331219299*w^5 - 408443308*w^4 - 151597163*w^3 - 408069089*w^2 + 227194757*w + 450126424");
    polCim3R.emplace_back("t+x+2*y+3*z+4*w");

    // Create cim3R ideal;
    Ideal<eltType> cim3R(polCim3R, 5 ,1000000);
    
    // Compute a reduced groebner basis;
    nbGen=cim3R.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        cim3R.printReducedGroebnerBasis("benchmark-cim3R.res", modulo);
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
    nbGen=cim3RF4(magma);

    cout << "cim3R: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


