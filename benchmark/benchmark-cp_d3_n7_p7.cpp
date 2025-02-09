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


int cp_d3_n7_p7F4(bool magma)
{
    cout << "#########################################################" << endl;
    cout << "#                      CP(3,7,7)                        #" << endl;
    cout << "#########################################################" << endl << endl;
    
    // Init element-prime tools
    eltType::setModulo(modulo);
    
    // Number of generator
    int nbGen;
    
    // Init monomial tools
    Monomial::initMonomial(7);
    
    // Create variable name array
    string * vars = new string[7];
    for(int i = 0; i < 7; i++)
      {
	vars[i]='x'+to_string(i+1);
      }
    Monomial::setVariable(vars);

    // Create polynomial array
    vector<Polynomial<eltType>> polCp_D3_N7_P7;
    
    // Fill the polynomial array
    polCp_D3_N7_P7.emplace_back("-7*x1^3+22*x1^2*x2+97*x1*x2^2+98*x2^3-55*x1^2*x3-73*x1*x2*x3-23*x2^2*x3+80*x1*x3^2+11*x2*x3^2-96*x3^3-94*x1^2*x4-4*x1*x2*x4+10*x2^2*x4-44*x1*x3*x4-49*x2*x3*x4+72*x3^2*x4-7*x1*x4^2+68*x2*x4^2+92*x3*x4^2-83*x4^3+87*x1^2*x5-83*x1*x2*x5-61*x2^2*x5+71*x1*x3*x5-47*x2*x3*x5-87*x3^2*x5-40*x1*x4*x5-10*x2*x4*x5-91*x3*x4*x5+98*x4^2*x5+75*x1*x5^2+95*x2*x5^2-28*x3*x5^2+37*x4*x5^2-90*x5^3-56*x1^2*x6-10*x1*x2*x6-8*x2^2*x6-17*x1*x3*x6+40*x2*x3*x6+47*x3^2*x6+42*x1*x4*x6+31*x2*x4*x6-88*x3*x4*x6-48*x4^2*x6-92*x1*x5*x6+x2*x5*x6+5*x3*x5*x6+5*x4*x5*x6-60*x5^2*x6+72*x1*x6^2-28*x2*x6^2-82*x3*x6^2+25*x4*x6^2+44*x5*x6^2-72*x6^3+62*x1*x2*x7-29*x2^2*x7-75*x1*x3*x7-81*x2*x3*x7-90*x3^2*x7-50*x1*x4*x7-51*x2*x4*x7-48*x3*x4*x7-19*x4^2*x7+6*x1*x5*x7+x2*x5*x7+13*x3*x5*x7+96*x4*x5*x7-34*x5^2*x7+37*x1*x6*x7+16*x2*x6*x7+71*x3*x6*x7+91*x4*x6*x7-2*x5*x6*x7-97*x6^2*x7+87*x1*x7^2-27*x2*x7^2+83*x3*x7^2+98*x4*x7^2-47*x5*x7^2+10*x6*x7^2+65*x7^3-62*x1^2-82*x1*x2+95*x2^2-10*x1*x3+91*x2*x3+43*x3^2+23*x1*x4+77*x2*x4+53*x3*x4+62*x4^2+74*x1*x5+55*x2*x5-10*x3*x5-17*x4*x5-13*x5^2-23*x1*x6+30*x2*x6+16*x3*x6+71*x5*x6+33*x6^2+44*x1*x7-15*x2*x7+9*x3*x7-64*x4*x7-39*x5*x7+7*x6*x7+12*x7^2+29*x1-59*x2-60*x3+64*x4-53*x5-89*x6-25*x7-96");
    polCp_D3_N7_P7.emplace_back("50*x1^3-60*x1^2*x2-60*x1*x2^2-10*x2^3-42*x1^2*x3+16*x1*x2*x3-65*x2^2*x3+69*x1*x3^2+36*x2*x3^2+57*x3^3+7*x1^2*x4+52*x1*x2*x4-85*x2^2*x4+80*x1*x3*x4+91*x2*x3*x4-49*x3^2*x4-35*x1*x4^2+60*x2*x4^2-29*x3*x4^2+90*x4^3-89*x1^2*x5-20*x1*x2*x5+54*x2^2*x5+28*x1*x3*x5-22*x2*x3*x5+31*x3^2*x5+97*x1*x4*x5-91*x2*x4*x5+5*x3*x4*x5+74*x4^2*x5-16*x1*x5^2-31*x2*x5^2+97*x3*x5^2+81*x4*x5^2-3*x5^3-70*x1^2*x6-4*x1*x2*x6-42*x1*x3*x6+51*x2*x3*x6+73*x3^2*x6+30*x1*x4*x6-47*x2*x4*x6-26*x3*x4*x6+27*x4^2*x6+59*x1*x5*x6+25*x2*x5*x6-67*x3*x5*x6+65*x4*x5*x6-56*x5^2*x6-33*x1*x6^2+65*x2*x6^2+37*x3*x6^2+5*x4*x6^2+42*x5*x6^2-51*x6^3+34*x1^2*x7-89*x1*x2*x7+18*x2^2*x7-33*x1*x3*x7-27*x2*x3*x7+95*x3^2*x7-64*x1*x4*x7-97*x2*x4*x7-51*x3*x4*x7+9*x4^2*x7-69*x1*x5*x7+31*x2*x5*x7+58*x3*x5*x7-12*x4*x5*x7-91*x5^2*x7+87*x1*x6*x7+88*x2*x6*x7+5*x3*x6*x7-63*x4*x6*x7+9*x5*x6*x7+16*x6^2*x7+40*x1*x7^2-6*x2*x7^2-57*x3*x7^2+36*x4*x7^2-27*x5*x7^2-44*x6*x7^2+49*x7^3-68*x1^2-77*x1*x2+52*x2^2+21*x1*x3+50*x2*x3+68*x3^2+89*x1*x4-2*x2*x4+88*x3*x4-91*x4^2-46*x1*x5-27*x2*x5+29*x3*x5+78*x4*x5-70*x5^2-34*x1*x6+10*x2*x6-36*x3*x6-5*x4*x6-21*x5*x6-85*x6^2+77*x1*x7+80*x2*x7+85*x3*x7-8*x4*x7-79*x5*x7-31*x6*x7-58*x7^2+x1-84*x2+80*x3+30*x4-22*x5+45*x6+49*x7+1");
    polCp_D3_N7_P7.emplace_back("-95*x1^3+86*x1^2*x2+62*x1*x2^2+7*x2^3-97*x1^2*x3+96*x1*x2*x3-61*x2^2*x3-95*x1*x3^2+46*x2*x3^2-52*x3^3-14*x1^2*x4-51*x1*x2*x4+45*x2^2*x4+61*x1*x3*x4-49*x2*x3*x4+37*x3^2*x4+57*x1*x4^2+38*x2*x4^2+63*x3*x4^2+83*x4^3+83*x1^2*x5+89*x1*x2*x5-50*x2^2*x5-2*x1*x3*x5-66*x2*x3*x5-59*x3^2*x5+28*x1*x4*x5-48*x2*x4*x5+21*x3*x4*x5-19*x4^2*x5-66*x1*x5^2-87*x2*x5^2+55*x3*x5^2-84*x4*x5^2+53*x5^3-96*x1^2*x6+14*x1*x2*x6-22*x2^2*x6+86*x1*x3*x6+18*x2*x3*x6-54*x3^2*x6+63*x1*x4*x6-12*x3*x4*x6-62*x4^2*x6-34*x1*x5*x6-63*x2*x5*x6-17*x3*x5*x6-83*x4*x5*x6-6*x5^2*x6-68*x1*x6^2+66*x2*x6^2-28*x3*x6^2+76*x4*x6^2+45*x5*x6^2+35*x6^3-8*x1^2*x7-79*x1*x2*x7+48*x2^2*x7+57*x1*x3*x7+16*x2*x3*x7-22*x3^2*x7+21*x1*x4*x7-87*x2*x4*x7-21*x3*x4*x7-12*x4^2*x7+72*x1*x5*x7-31*x2*x5*x7+83*x3*x5*x7+96*x4*x5*x7+52*x5^2*x7-15*x1*x6*x7+88*x2*x6*x7-5*x3*x6*x7+95*x4*x6*x7-86*x5*x6*x7-31*x6^2*x7+17*x1*x7^2-39*x2*x7^2+61*x4*x7^2-21*x5*x7^2+17*x6*x7^2+44*x7^3-54*x1^2-58*x1*x2-82*x2^2-35*x1*x3-22*x2*x3+17*x3^2-71*x1*x4-13*x2*x4-74*x3*x4+48*x4^2-40*x1*x5+90*x2*x5-65*x3*x5+6*x4*x5-40*x5^2-32*x1*x6+74*x2*x6+27*x3*x6-14*x4*x6-88*x5*x6-17*x6^2+87*x1*x7-20*x2*x7+26*x3*x7+83*x4*x7+28*x5*x7+31*x6*x7-16*x7^2-60*x1+31*x2+13*x3-83*x4+92*x5-52*x6-37*x7-59");
    polCp_D3_N7_P7.emplace_back("-47*x1^3+58*x1^2*x2+93*x1*x2^2+18*x2^3-56*x1^2*x3+74*x1*x2*x3+75*x2^2*x3-x1*x3^2-91*x2*x3^2+8*x3^3-45*x1^2*x4+95*x1*x2*x4+10*x2^2*x4+49*x1*x3*x4-53*x2*x3*x4+43*x3^2*x4+46*x1*x4^2-72*x2*x4^2-2*x3*x4^2+42*x4^3-15*x1^2*x5+65*x1*x2*x5+18*x2^2*x5-27*x1*x3*x5+37*x2*x3*x5+58*x3^2*x5-36*x1*x4*x5+2*x2*x4*x5-45*x3*x4*x5+86*x4^2*x5-34*x1*x5^2-83*x2*x5^2-77*x3*x5^2-4*x4*x5^2-50*x5^3-67*x1^2*x6-10*x1*x2*x6+60*x2^2*x6+88*x1*x3*x6+76*x2*x3*x6-84*x3^2*x6-76*x1*x4*x6+4*x2*x4*x6-85*x3*x4*x6+7*x4^2*x6+52*x1*x5*x6-45*x3*x5*x6+5*x4*x5*x6-75*x5^2*x6+11*x1*x6^2-46*x2*x6^2+67*x3*x6^2+4*x4*x6^2+14*x5*x6^2+53*x6^3-62*x1^2*x7+11*x1*x2*x7-90*x2^2*x7-63*x1*x3*x7-44*x2*x3*x7+78*x3^2*x7-22*x1*x4*x7-84*x2*x4*x7+58*x3*x4*x7-11*x4^2*x7+x1*x5*x7+84*x2*x5*x7-34*x3*x5*x7-40*x4*x5*x7+11*x5^2*x7-11*x1*x6*x7+63*x2*x6*x7+95*x3*x6*x7-89*x4*x6*x7-56*x5*x6*x7-68*x6^2*x7+94*x1*x7^2-43*x2*x7^2+45*x3*x7^2-81*x4*x7^2+80*x5*x7^2-5*x6*x7^2+24*x7^3+70*x1^2-56*x1*x2+11*x2^2+61*x1*x3+67*x2*x3-12*x3^2-55*x1*x4-60*x2*x4-62*x3*x4+2*x4^2+87*x1*x5+77*x2*x5+94*x3*x5+5*x4*x5+80*x5^2+88*x1*x6+27*x2*x6-38*x3*x6-54*x4*x6+8*x5*x6+27*x6^2+6*x1*x7+72*x2*x7-99*x3*x7+22*x4*x7-91*x5*x7+33*x6*x7-18*x7^2+21*x1-64*x2-80*x3-38*x4-48*x5-87*x6-8*x7+37");
    polCp_D3_N7_P7.emplace_back("-96*x1^3-82*x1^2*x2-86*x1*x2^2+67*x2^3+50*x1^2*x3-80*x1*x2*x3+39*x2^2*x3+5*x1*x3^2+81*x2*x3^2+2*x3^3+88*x1^2*x4-96*x1*x2*x4-27*x2^2*x4+65*x1*x3*x4+39*x2*x3*x4+7*x3^2*x4-21*x1*x4^2-10*x2*x4^2+44*x3*x4^2-42*x4^3-73*x1^2*x5-72*x1*x2*x5+11*x2^2*x5-45*x1*x3*x5+76*x2*x3*x5+24*x3^2*x5-80*x1*x4*x5+61*x2*x4*x5+41*x3*x4*x5-97*x4^2*x5-28*x1*x5^2-4*x2*x5^2-64*x3*x5^2+89*x4*x5^2-52*x5^3-88*x1^2*x6-16*x1*x2*x6+72*x2^2*x6-59*x1*x3*x6+57*x2*x3*x6-91*x3^2*x6+98*x1*x4*x6-57*x2*x4*x6-65*x3*x4*x6-78*x4^2*x6-47*x1*x5*x6-76*x2*x5*x6-62*x3*x5*x6-3*x4*x5*x6+51*x5^2*x6+31*x1*x6^2+64*x2*x6^2+11*x3*x6^2-89*x4*x6^2-27*x5*x6^2-39*x6^3-52*x1^2*x7-11*x1*x2*x7-51*x2^2*x7+83*x1*x3*x7+78*x2*x3*x7-16*x3^2*x7+77*x1*x4*x7+26*x2*x4*x7-20*x3*x4*x7+17*x4^2*x7+34*x1*x5*x7+92*x2*x5*x7-86*x3*x5*x7+80*x4*x5*x7-75*x5^2*x7-81*x1*x6*x7-21*x2*x6*x7+78*x3*x6*x7+58*x4*x6*x7-38*x5*x6*x7+44*x6^2*x7-15*x1*x7^2+86*x2*x7^2-63*x3*x7^2-61*x4*x7^2-77*x5*x7^2+75*x6*x7^2-80*x7^3+3*x1^2+94*x1*x2-33*x2^2-3*x1*x3+49*x2*x3-98*x3^2-60*x1*x4-8*x2*x4-55*x3*x4-43*x4^2-4*x1*x5-2*x2*x5+67*x3*x5-37*x4*x5+23*x5^2+98*x1*x6-81*x2*x6+38*x3*x6+32*x4*x6-14*x5*x6-49*x6^2+42*x1*x7+71*x2*x7+89*x3*x7+83*x4*x7+24*x5*x7+14*x6*x7-17*x7^2-37*x1-80*x2+30*x3-9*x4+4*x5+61*x6-90*x7+13");
    polCp_D3_N7_P7.emplace_back("68*x1^3-84*x1^2*x2-95*x1*x2^2-25*x2^3-81*x1^2*x3-22*x1*x2*x3-26*x2^2*x3-64*x1*x3^2-67*x2*x3^2-x3^3-60*x1^2*x4+57*x1*x2*x4-7*x2^2*x4+70*x1*x3*x4+88*x2*x3*x4+23*x3^2*x4+x1*x4^2+7*x2*x4^2-14*x3*x4^2-48*x4^3+68*x1^2*x5+18*x1*x2*x5-93*x2^2*x5+97*x1*x3*x5-5*x2*x3*x5-83*x3^2*x5+11*x1*x4*x5-42*x2*x4*x5-92*x3*x4*x5+51*x4^2*x5-63*x1*x5^2+55*x2*x5^2-57*x3*x5^2-5*x4*x5^2+3*x5^3+11*x1^2*x6+51*x1*x2*x6-36*x2^2*x6+21*x1*x3*x6-36*x2*x3*x6-26*x3^2*x6+72*x1*x4*x6-73*x2*x4*x6-90*x3*x4*x6-56*x4^2*x6-63*x1*x5*x6-49*x2*x5*x6+93*x3*x5*x6+7*x4*x5*x6+6*x5^2*x6-56*x1*x6^2+17*x2*x6^2-28*x3*x6^2-31*x4*x6^2+93*x5*x6^2-40*x6^3-29*x1^2*x7+49*x1*x2*x7+53*x2^2*x7-27*x1*x3*x7-x2*x3*x7+35*x3^2*x7-56*x1*x4*x7+69*x2*x4*x7-32*x3*x4*x7-20*x4^2*x7+94*x1*x5*x7+86*x2*x5*x7+82*x3*x5*x7-15*x4*x5*x7+60*x5^2*x7+73*x1*x6*x7+31*x2*x6*x7-90*x3*x6*x7+57*x4*x6*x7-80*x5*x6*x7-82*x6^2*x7+80*x1*x7^2+36*x2*x7^2-45*x3*x7^2-80*x4*x7^2-40*x5*x7^2-69*x6*x7^2-37*x7^3-25*x1^2+34*x1*x2-78*x2^2-9*x1*x3+31*x2*x3+34*x3^2+53*x1*x4+24*x2*x4+15*x3*x4-77*x4^2+89*x1*x5-23*x2*x5-85*x3*x5+62*x4*x5+17*x5^2-x1*x6+75*x2*x6-27*x3*x6+98*x4*x6+88*x5*x6-82*x6^2+69*x1*x7+14*x2*x7+20*x3*x7-51*x4*x7-10*x5*x7+3*x6*x7+42*x7^2+76*x1-55*x2-49*x3+39*x4+67*x5+69*x6-51*x7-80");
    polCp_D3_N7_P7.emplace_back("-72*x1^3+36*x1^2*x2+50*x1*x2^2-37*x2^3+71*x1^2*x3-19*x1*x2*x3-34*x2^2*x3+83*x1*x3^2-43*x2*x3^2-27*x3^3+54*x1^2*x4-18*x1*x2*x4-12*x2^2*x4-15*x1*x3*x4-34*x2*x3*x4+78*x3^2*x4+60*x1*x4^2+14*x2*x4^2+13*x3*x4^2+84*x1^2*x5-45*x1*x2*x5+9*x2^2*x5-80*x1*x3*x5-20*x2*x3*x5+25*x3^2*x5+80*x1*x4*x5-27*x2*x4*x5+88*x3*x4*x5-40*x4^2*x5+82*x1*x5^2-42*x2*x5^2+50*x3*x5^2-57*x4*x5^2-88*x5^3-53*x1^2*x6-30*x1*x2*x6+63*x2^2*x6-40*x1*x3*x6+78*x2*x3*x6+30*x3^2*x6-33*x1*x4*x6-27*x2*x4*x6-96*x3*x4*x6+48*x4^2*x6+83*x1*x5*x6+30*x2*x5*x6-55*x3*x5*x6-75*x4*x5*x6-17*x5^2*x6+97*x1*x6^2+13*x2*x6^2+90*x3*x6^2-84*x4*x6^2-75*x5*x6^2-96*x6^3+65*x1^2*x7+73*x1*x2*x7+16*x2^2*x7+18*x1*x3*x7-98*x2*x3*x7+17*x3^2*x7+67*x1*x4*x7-99*x2*x4*x7-67*x3*x4*x7-12*x4^2*x7+94*x1*x5*x7-21*x2*x5*x7-64*x3*x5*x7-71*x4*x5*x7-32*x5^2*x7-76*x1*x6*x7+14*x2*x6*x7-80*x3*x6*x7-47*x4*x6*x7+82*x5*x6*x7-6*x6^2*x7-32*x1*x7^2-77*x2*x7^2+36*x3*x7^2+21*x4*x7^2-19*x5*x7^2+x6*x7^2-58*x7^3-9*x1^2-37*x1*x2+7*x2^2+64*x1*x3-86*x2*x3-65*x3^2-14*x1*x4-81*x2*x4-85*x3*x4-88*x4^2+77*x1*x5-27*x2*x5-75*x3*x5+21*x4*x5-39*x5^2+19*x1*x6-9*x2*x6+91*x3*x6+63*x4*x6-93*x5*x6+x6^2+79*x1*x7-35*x2*x7-86*x3*x7-66*x4*x7+3*x5*x7+62*x6*x7-66*x7^2-4*x1+29*x2+97*x3-24*x4+65*x5-25*x6+56*x7+3");

    // Create cp_d3_n7_p7 ideal;
    Ideal<eltType> cp_d3_n7_p7(polCp_D3_N7_P7, 7, 20000000);
    
    // Compute a reduced groebner basis;
    nbGen=cp_d3_n7_p7.f4();
    
    // Print the reduced groebner basis into a file
    if(magma)
    {
        cp_d3_n7_p7.printReducedGroebnerBasis("benchmark-cp_d3_n7_p7.res", modulo);
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
    nbGen=cp_d3_n7_p7F4(magma);

    cout << "cp(3,7,7): " << chrono::duration_cast<millisecs_t>(chrono::steady_clock::now()-start).count() << " ms                   (" << nbGen << " generators)" << endl << endl;
   
    return 0;
}


