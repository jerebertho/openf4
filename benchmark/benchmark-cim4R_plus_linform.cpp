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


int cim4RF4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                         CIM4R                         #" << endl;
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
    vector<Polynomial<eltType>> polCim4R;
    
    // Fill the polynomial array
    polCim4R.emplace_back("219201386*x^8 - 459910334*x^7 + 33367490*x^6 - 285538596*x^5 - 377848199*x^4 - 68051022*x^3 + 236915708*x^2 - 156254838*x + 197223747*y^8 - 20168465*y^7 + 153650995*y^6 - 174872719*y^5 + 183825463*y^4 - 233215643*y^3 - 224255198*y^2 + 440140637*y - 219201386*z^8 - 200709750*z^7 - 142697255*z^6 - 49118975*z^5 + 137665636*z^4 - 81056313*z^3 - 482453417*z^2 + 128852289*z - 197223747*w^8 + 108479711*w^7 + 168924457*w^6 - 219488832*w^5 + 205972622*w^4 - 257300335*w^3 + 16847106*w^2 - 110327721*w + 196212944");
    polCim4R.emplace_back("-197223747*x^8 + 217282826*x^7 - 181737589*x^6 + 23794883*x^5 - 226538785*x^4 - 320102451*x^3 - 509140023*x^2 - 8466177*x + 137395811*y^16 - 482123228*y^15 - 9329438*y^14 + 210428282*y^13 + 365050461*y^12 + 455369960*y^11 - 88214720*y^10 + 215534968*y^9 - 46064988*y^8 + 332150528*y^7 - 471630642*y^6 + 393134240*y^5 - 111793740*y^4 - 462018460*y^3 + 503833225*y^2 + 360490819*y + 197223747*z^8 - 520145082*z^7 + 459283850*z^6 + 451485660*z^5 + 336032481*z^4 - 523992741*z^3 - 422719311*z^2 - 503355896*z - 137395811*w^16 + 238020848*w^15 + 107392677*w^14 - 269616668*w^13 - 191219039*w^12 - 238788872*w^11 + 199197217*w^10 - 304300052*w^9 + 475238625*w^8 + 304436464*w^7 + 332383866*w^6 - 484312229*w^5 - 24594998*w^4 - 22127869*w^3 + 509228262*w^2 - 300733812*w - 40290953");
    polCim4R.emplace_back("-137395811*x^16 + 426752524*x^15 - 441761473*x^14 + 532857448*x^13 - 366888483*x^12 + 405780818*x^11 + 398379959*x^10 - 404815674*x^9 + 490844352*x^8 - 332624472*x^7 - 508390739*x^6 - 490382320*x^5 + 397686558*x^4 - 29389038*x^3 + 267425701*x^2 - 472737266*x + 81198893*y^24 + 287402795*y^23 + 225162610*y^22 + 529947542*y^21 - 108621467*y^20 + 152639741*y^19 - 297712296*y^18 + 306179429*y^17 + 440644089*y^16 + 335921699*y^15 + 197694272*y^14 - 280226100*y^13 + 392347812*y^12 + 334563197*y^11 + 250461584*y^10 - 402475424*y^9 - 137904221*y^8 + 535844711*y^7 + 186504127*y^6 - 385410022*y^5 - 523083921*y^4 - 462068873*y^3 + 199981049*y^2 + 463869835*y + 243172066");
    polCim4R.emplace_back("-137395811*z^16 - 242668934*z^15 - 392861104*z^14 - 515047616*z^13 + 7959987*z^12 + 123008503*z^11 + 52067084*z^10 - 79372358*z^9 - 285231505*z^8 - 485449285*z^7 + 160510508*z^6 + 175854858*z^5 + 391293704*z^4 + 367370727*z^3 - 329805037*z^2 + 320631281*z + 81198893*w^24 - 408945137*w^23 + 470083978*w^22 - 274353201*w^21 + 82384644*w^20 + 90200366*w^19 + 287950097*w^18 + 115810603*w^17 + 257657773*w^16 - 311450646*w^15 + 444586590*w^14 + 67193455*w^13 + 530458300*w^12 - 257024998*w^11 + 12984029*w^10 + 474790861*w^9 - 506168320*w^8 - 488167963*w^7 - 294072668*w^6 - 441297916*w^5 - 448824097*w^4 - 285290981*w^3 - 448446964*w^2 + 256451518*w - 394953658");
    polCim4R.emplace_back("t+x+2*y+3*z+4*w");

    // Create cim4R ideal;
    Ideal<eltType> cim4R(polCim4R, 5 ,1000000);
    
    // Compute a reduced groebner basis;
    nbGen=cim4R.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        cim4R.printReducedGroebnerBasis("benchmark-cim4R.res", modulo);
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
    nbGen=cim4RF4(magma);

    cout << "cim4R: " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;

    return 0;
}


