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


int cp_d4_n6_p6F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                      CP(4,6,6)                        #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(6);
    
    // Create variable name array
    string * vars = new string[6];
    for(int i = 0; i < 6; i++)
      {
	vars[i]='x'+to_string(i+1);
      }
    Monomial::setVariable(vars);

    // Create polynomial array
    vector<Polynomial<eltType>> polCp_D4_N6_P6;
    
    // Fill the polynomial array
    polCp_D4_N6_P6.emplace_back("1073741799*x1^3+66*x1^2*x2+1073741703*x1*x2^2+6*x2^3+1073741662*x1^2*x3+194*x1*x2*x3+74*x2^2*x3+124*x1*x3^2+44*x2*x3^2+68*x3^3+1073741545*x1^2*x4+1073741681*x1*x2*x4+72*x2^2*x4+1073741663*x1*x3*x4+29*x2*x3*x4+1073741817*x3^2*x4+1073741793*x1*x4^2+1073741766*x2*x4^2+95*x3*x4^2+1073741731*x4^3+261*x1^2*x5+1073741819*x1*x2*x5+37*x2^2*x5+160*x1*x3*x5+98*x2*x3*x5+31*x3^2*x5+1073741677*x1*x4*x5+1073741819*x2*x4*x5+x3*x4*x5+72*x4^2*x5+1073741747*x1*x5^2+11*x2*x5^2+1073741799*x3*x5^2+1073741737*x4*x5^2+53*x5^3+1073741659*x1^2*x6+1073741661*x1*x2*x6+1073741804*x2^2*x6+1073741739*x1*x3*x6+1073741804*x2*x3*x6+1073741776*x3^2*x6+1073741807*x1*x4*x6+1073741798*x2*x4*x6+x3*x4*x6+1073741740*x4^2*x6+84*x1*x5*x6+1073741778*x2*x5*x6+16*x3*x5*x6+43*x4*x5*x6+1073741799*x5^2*x6+46*x1*x6^2+40*x2*x6^2+1073741800*x3*x6^2+1073741736*x4*x6^2+13*x5*x6^2+71*x6^3+1073741807*x1*x2+87*x2^2+142*x1*x3+10*x2*x3+77*x3^2+1073741813*x1*x4+95*x2*x4+55*x3*x4+47*x4^2+1073741727*x1*x5+1073741780*x2*x5+30*x3*x5+92*x4*x5+5*x5^2+150*x1*x6+1073741746*x2*x6+1073741812*x3*x6+1073741739*x4*x6+1073741817*x5*x6+16*x6^2+1073741643*x1+91*x2+1073741768*x3+1073741779*x4+1073741745*x5+83*x6+9");
    polCp_D4_N6_P6.emplace_back("22*x1^3+1073741703*x1^2*x2+18*x1*x2^2+1073741587*x2^3+97*x1^2*x3+148*x1*x2*x3+1073741578*x2^2*x3+44*x1*x3^2+74*x2*x3^2+1073741825*x3^3+1073741754*x1^2*x4+144*x1*x2*x4+294*x2^2*x4+29*x1*x3*x4+10*x2*x3*x4+71*x3^2*x4+1073741766*x1*x4^2+182*x2*x4^2+1073741755*x3*x4^2+50*x4^3+1073741823*x1^2*x5+74*x1*x2*x5+1073741683*x2^2*x5+98*x1*x3*x5+192*x2*x3*x5+1073741780*x3^2*x5+1073741819*x1*x4*x5+1073741730*x3*x4*x5+1073741767*x4^2*x5+11*x1*x5^2+128*x2*x5^2+7*x3*x5^2+1073741738*x4*x5^2+52*x5^3+1073741744*x1^2*x6+1073741781*x1*x2*x6+1073741770*x2^2*x6+1073741804*x1*x3*x6+1073741793*x2*x3*x6+1073741788*x3^2*x6+1073741798*x1*x4*x6+196*x2*x4*x6+33*x3*x4*x6+1073741785*x4^2*x6+1073741778*x1*x5*x6+1073741647*x2*x5*x6+1073741738*x3*x5*x6+1073741757*x4*x5*x6+1073741807*x5^2*x6+40*x1*x6^2+1073741759*x2*x6^2+12*x3*x6^2+1073741759*x4*x6^2+1073741738*x5*x6^2+80*x6^3+1073741817*x1^2+174*x1*x2+186*x2^2+10*x1*x3+50*x2*x3+1073741774*x3^2+95*x1*x4+1073741699*x2*x4+10*x3*x4+7*x4^2+1073741780*x1*x5+1073741707*x2*x5+65*x3*x5+34*x4*x5+1073741823*x5^2+1073741746*x1*x6+1073741801*x2*x6+1073741802*x3*x6+1073741767*x4*x6+1073741750*x5*x6+28*x6^2+91*x1+88*x2+1073741731*x3+16*x4+69*x5+1073741785*x6+1073741794");
    polCp_D4_N6_P6.emplace_back("1073741772*x1^3+97*x1^2*x2+74*x1*x2^2+1073741744*x2^3+124*x1^2*x3+88*x1*x2*x3+74*x2^2*x3+204*x1*x3^2+1073741821*x2*x3^2+84*x3^3+1073741745*x1^2*x4+29*x1*x2*x4+5*x2^2*x4+1073741807*x1*x3*x4+142*x2*x3*x4+1073741722*x3^2*x4+95*x1*x4^2+1073741755*x2*x4^2+178*x3*x4^2+x4^3+80*x1^2*x5+98*x1*x2*x5+96*x2^2*x5+62*x1*x3*x5+1073741733*x2*x3*x5+291*x3^2*x5+x1*x4*x5+1073741730*x2*x4*x5+1073741795*x3*x4*x5+1073741817*x4^2*x5+1073741799*x1*x5^2+7*x2*x5^2+1073741735*x3*x5^2+54*x4*x5^2+1073741805*x5^3+1073741783*x1^2*x6+1073741804*x1*x2*x6+1073741810*x2^2*x6+1073741725*x1*x3*x6+1073741749*x2*x3*x6+90*x3^2*x6+x1*x4*x6+33*x2*x4*x6+118*x3*x4*x6+1073741762*x4^2*x6+16*x1*x5*x6+1073741738*x2*x5*x6+1073741761*x3*x5*x6+51*x5^2*x6+1073741800*x1*x6^2+12*x2*x6^2+1073741759*x3*x6^2+52*x4*x6^2+50*x5*x6^2+1073741780*x6^3+71*x1^2+10*x1*x2+25*x2^2+154*x1*x3+1073741721*x2*x3+1073741635*x3^2+55*x1*x4+10*x2*x4+1073741689*x3*x4+1073741742*x4^2+30*x1*x5+65*x2*x5+174*x3*x5+18*x4*x5+1073741800*x5^2+1073741812*x1*x6+1073741802*x2*x6+80*x3*x6+36*x4*x6+60*x5*x6+1073741730*x6^2+1073741768*x1+1073741731*x2+154*x3+91*x4+1073741736*x5+1073741825*x6+1073741796");
    polCp_D4_N6_P6.emplace_back("1073741733*x1^3+1073741754*x1^2*x2+72*x1*x2^2+98*x2^3+1073741745*x1^2*x3+29*x1*x2*x3+5*x2^2*x3+1073741817*x1*x3^2+71*x2*x3^2+1073741792*x3^3+1073741793*x1^2*x4+1073741705*x1*x2*x4+182*x2^2*x4+190*x1*x3*x4+1073741683*x2*x3*x4+178*x3^2*x4+1073741539*x1*x4^2+150*x2*x4^2+3*x3*x4^2+100*x4^3+1073741752*x1^2*x5+1073741819*x1*x2*x5+x1*x3*x5+1073741730*x2*x3*x5+1073741811*x3^2*x5+144*x1*x4*x5+1073741707*x2*x4*x5+1073741807*x3*x4*x5+93*x4^2*x5+1073741737*x1*x5^2+1073741738*x2*x5^2+54*x3*x5^2+176*x4*x5^2+1073741778*x5^3+1073741817*x1^2*x6+1073741798*x1*x2*x6+98*x2^2*x6+x1*x3*x6+33*x2*x3*x6+59*x3^2*x6+1073741653*x1*x4*x6+1073741743*x2*x4*x6+1073741697*x3*x4*x6+1073741746*x4^2*x6+43*x1*x5*x6+1073741757*x2*x5*x6+20*x4*x5*x6+31*x5^2*x6+1073741736*x1*x6^2+1073741759*x2*x6^2+52*x3*x6^2+160*x4*x6^2+95*x5*x6^2+5*x6^3+1073741820*x1^2+95*x1*x2+1073741763*x2^2+55*x1*x3+10*x2*x3+1073741758*x3^2+94*x1*x4+14*x2*x4+1073741657*x3*x4+195*x4^2+92*x1*x5+34*x2*x5+18*x3*x5+1073741815*x4*x5+73*x5^2+1073741739*x1*x6+1073741767*x2*x6+36*x3*x6+1073741659*x4*x6+68*x5*x6+1073741801*x6^2+1073741779*x1+16*x2+91*x3+114*x4+1073741798*x5+1073741776*x6+88");
    polCp_D4_N6_P6.emplace_back("87*x1^3+1073741823*x1^2*x2+37*x1*x2^2+1073741779*x2^3+80*x1^2*x3+98*x1*x2*x3+96*x2^2*x3+31*x1*x3^2+1073741780*x2*x3^2+97*x3^3+1073741752*x1^2*x4+1073741819*x1*x2*x4+x1*x3*x4+1073741730*x2*x3*x4+1073741811*x3^2*x4+72*x1*x4^2+1073741767*x2*x4^2+1073741817*x3*x4^2+31*x4^3+1073741747*x1^2*x5+22*x1*x2*x5+128*x2^2*x5+1073741771*x1*x3*x5+14*x2*x3*x5+1073741735*x3^2*x5+1073741647*x1*x4*x5+1073741649*x2*x4*x5+108*x3*x4*x5+176*x4^2*x5+159*x1*x5^2+156*x2*x5^2+1073741761*x3*x5^2+1073741680*x4*x5^2+388*x5^3+42*x1^2*x6+1073741778*x1*x2*x6+1073741737*x2^2*x6+16*x1*x3*x6+1073741738*x2*x3*x6+1073741794*x3^2*x6+43*x1*x4*x6+1073741757*x2*x4*x6+10*x4^2*x6+1073741771*x1*x5*x6+1073741787*x2*x5*x6+102*x3*x5*x6+62*x4*x5*x6+1073741626*x5^2*x6+13*x1*x6^2+1073741738*x2*x6^2+50*x3*x6^2+95*x4*x6^2+58*x5*x6^2+1073741791*x6^3+1073741777*x1^2+1073741780*x1*x2+1073741767*x2^2+30*x1*x3+65*x2*x3+87*x3^2+92*x1*x4+34*x2*x4+18*x3*x4+1073741821*x4^2+10*x1*x5+1073741819*x2*x5+1073741773*x3*x5+146*x4*x5+174*x5^2+1073741817*x1*x6+1073741750*x2*x6+60*x3*x6+68*x4*x6+74*x5*x6+1073741770*x6^2+1073741745*x1+69*x2+1073741736*x3+1073741798*x4+10*x5+85*x6+80");
    polCp_D4_N6_P6.emplace_back("1073741771*x1^3+1073741744*x1^2*x2+1073741804*x1*x2^2+1073741808*x2^3+1073741783*x1^2*x3+1073741804*x1*x2*x3+1073741810*x2^2*x3+1073741776*x1*x3^2+1073741788*x2*x3^2+30*x3^3+1073741817*x1^2*x4+1073741798*x1*x2*x4+98*x2^2*x4+x1*x3*x4+33*x2*x3*x4+59*x3^2*x4+1073741740*x1*x4^2+1073741785*x2*x4^2+1073741762*x3*x4^2+1073741800*x4^3+42*x1^2*x5+1073741778*x1*x2*x5+1073741737*x2^2*x5+16*x1*x3*x5+1073741738*x2*x3*x5+1073741794*x3^2*x5+43*x1*x4*x5+1073741757*x2*x4*x5+10*x4^2*x5+1073741799*x1*x5^2+1073741807*x2*x5^2+51*x3*x5^2+31*x4*x5^2+1073741760*x5^3+46*x1^2*x6+80*x1*x2*x6+1073741759*x2^2*x6+1073741773*x1*x3*x6+24*x2*x3*x6+1073741759*x3^2*x6+1073741645*x1*x4*x6+1073741691*x2*x4*x6+104*x3*x4*x6+160*x4^2*x6+26*x1*x5*x6+1073741649*x2*x5*x6+100*x3*x5*x6+190*x4*x5*x6+58*x5^2*x6+213*x1*x6^2+240*x2*x6^2+1073741686*x3*x6^2+15*x4*x6^2+1073741719*x5*x6^2+360*x6^3+75*x1^2+1073741746*x1*x2+1073741814*x2^2+1073741812*x1*x3+1073741802*x2*x3+40*x3^2+1073741739*x1*x4+1073741767*x2*x4+36*x3*x4+1073741743*x4^2+1073741817*x1*x5+1073741750*x2*x5+60*x3*x5+68*x4*x5+37*x5^2+32*x1*x6+56*x2*x6+1073741633*x3*x6+1073741775*x4*x6+1073741713*x5*x6+222*x6^2+83*x1+1073741785*x2+1073741825*x3+1073741776*x4+85*x5+54*x6+9");

    // Create cp_d4_n6_p6 ideal;
    Ideal<eltType> cp_d4_n6_p6(polCp_D4_N6_P6, 6, 20000000);
    
    // Compute a reduced groebner basis;
    nbGen=cp_d4_n6_p6.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        cp_d4_n6_p6.printReducedGroebnerBasis("benchmark-cp_d4_n6_p6.res", modulo);
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
    nbGen=cp_d4_n6_p6F4(magma);

    cout << "cp(4,6,6): " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;
   
    return 0;
}


