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


int cim3_5RF4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         CIM3_5R                       #" << endl;
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
    vector<Polynomial<eltType>> polCim3_5R;
    
    // Fill the polynomial array
    polCim3_5R.emplace_back("219201386*x^8 + 81697143*x^7 - 486640724*x^6 - 172417078*x^5 - 299572653*x^4 - 440028745*x^3 - 146217024*x^2 - 333377164*x - 20111*y^4 + 256389651*y^3 - 176186055*y^2 + 449572913*y - 219201386*z^8 + 381996515*z^7 + 89080350*z^6 + 351316560*z^5 + 286814791*z^4 + 156526964*z^3 + 493248253*z^2 + 259719573*z + 20111*w^4 + 477655188*w^3 + 332878824*w^2 - 141845088*w - 353054872");
    polCim3_5R.emplace_back("-197223747*x^8 + 261627352*x^7 - 317360536*x^6 - 164588875*x^5 - 412018712*x^4 - 244013652*x^3 + 441328817*x^2 + 180386356*x + 515281093*y^8 + 239586225*y^7 + 272284570*y^6 + 68394118*y^5 - 124724710*y^4 - 119709512*y^3 - 301428657*y^2 + 146308135*y + 197223747*z^8 - 479225455*z^7 + 529447532*z^6 - 29276380*z^5 - 481702957*z^4 + 193789563*z^3 - 496494025*z^2 - 120933485*z - 515281093*w^8 - 430822514*w^7 + 415765871*w^6 + 139599093*w^5 - 157295429*w^4 - 151925913*w^3 - 426313531*w^2 + 69291238*w + 200934145");
    polCim3_5R.emplace_back("-137395811*x^16 + 262564210*x^15 + 7931417*x^14 - 279925877*x^13 + 385653133*x^12 - 141280232*x^11 - 426548974*x^10 - 301737027*x^9 + 260498629*x^8 + 231786892*x^7 - 177403202*x^6 - 181964954*x^5 - 357982989*x^4 + 48877749*x^3 - 181334478*x^2 - 372670493*x - 136055684*y^12 + 414996191*y^11 - 150322389*y^10 - 336140809*y^9 + 247829001*y^8 - 355584738*y^7 - 468926186*y^6 + 410105986*y^5 + 144142871*y^4 + 213754270*y^3 - 88430069*y^2 - 295675439*y - 430910171");
    polCim3_5R.emplace_back("-137395811*z^16 - 209642085*z^15 - 175090601*z^14 + 80999313*z^13 - 324947844*z^12 + 129122993*z^11 - 413796337*z^10 + 89129155*z^9 + 501892283*z^8 + 381783820*z^7 + 274690460*z^6 + 194463995*z^5 + 244241111*z^4 - 469408155*z^3 + 485135825*z^2 - 454104978*z - 136055684*w^12 + 163245475*w^11 + 59138502*w^10 + 385670862*w^9 + 53684269*w^8 - 243911227*w^7 + 344199548*w^6 - 169513668*w^5 + 109100299*w^4 - 165726602*w^3 + 529929221*w^2 + 39752519*w + 326877609");
    polCim3_5R.emplace_back("t+x+2*y+3*z+4*w");

    // Create cim3_5R ideal;
    Ideal<eltType> cim3_5R(polCim3_5R, 5 ,1000000);
    
    // Compute a reduced groebner basis;
    nbGen=cim3_5R.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        cim3_5R.printReducedGroebnerBasis("benchmark-cim3_5R.res", modulo);
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
    nbGen=cim3_5RF4(magma);

    cout << "cim3_5R: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


