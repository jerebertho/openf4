/* 
 * Copyright (C) 2010 Antoine Joux and Vanessa Vitse 

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

  /**
  * \file polynomial.h
  * \brief Declaration of class Polynomial.
  * \author Vanessa VITSE, Antoine JOUX, Titouan COLADON
  */

#ifndef F4_POLYNOMIAL_H
#define F4_POLYNOMIAL_H

#include <forward_list>
#include "term.h"

/** \namespace F4 
 * Group all the required tools used by the F4 algorithm.
 */
namespace F4
{
    /**
     * \class Polynomial
     * Represent a polynomial.
     */
    template <typename Element>
    class Polynomial
    {
        public:
        
            // Constructor
            
            /**
             * \brief Constructor.
             */
            Polynomial();
            
            /**
             * \brief Constructor.
             * \param s: String representing the polynomial.
             */
            Polynomial(std::string const s);
            
            /**
             * \brief Copy constructor.
             * \param polynomial: Polynomial to copy.
             */
            Polynomial(Polynomial const & polynomial);
            
            /**
             * \brief Move constructor.
             * \param polynomial: Polynomial to move.
             */
            Polynomial(Polynomial && polynomial);
            
            
            // Destructor
            
            /**
             * \brief Destructor.
             */
            ~Polynomial();
            
            
            // Miscellaneous
            
            /**
             * \brief Print the polynomial.
             */
            void printPolynomial (std::ostream & stream = std::cout) const;
            
            /**
             * \brief Get the number of terms of this.
             * \return Number of terms of this.
             */
            int getNbTerm() const;
             
            /**
            * \brief Get the leading term of this.
            * \pre _polynomial is not empty.
            * \return Leading term of this.
            */
            const Term<Element> & getLT() const; 
            
            /**
            * \brief Get the number of the leading monomial of this.
            * \pre _polynomial is not empty.
            * \return Number of the leading monomial of this.
            */
            int getLM() const;
            
            /**
            * \brief Get the leading coefficient of this.
            * \pre _polynomial is not empty.
            * \return Leading coefficient of this.
            */
            Element getLC() const;
            
            /**
            * \brief Get the coefficient of the term of monomial numMon.
            * \return Coefficient of the term of monomial numMon.
            */
            Element getCoefficient(int numMon) const;
            
            /**
            * \brief Delete the leading term of this.
            */
            void deleteLT(); 
            
            /**
             * \brief Reset this.
             */
            void reset(); 
            
            /**
             * \brief Normalize this.
             */
            //void normalize();
            
            
            // Internal operators
            
            /**
             * \brief Overload the operator =.
             * \param polynomial: Polynomial to copy.
             * \return Reference on this.
             */
            Polynomial & operator=(Polynomial const & polynomial);
            
            /**
             * \brief Overload the move operator =.
             * \param polynomial: Polynomial to move.
             * \return Reference on this.
             */
            Polynomial & operator=(Polynomial && polynomial);
            
            /**
             * \brief Overload the operator *= to multiply this with a monomial.
             * \param mon: Monomial.
             * \return Reference on this.
             */
            Polynomial & operator*=(Monomial const & monomial);
            
            /**
             * \brief Overload the operator *= to multiply this with a monomial.
             * \param numMon: number of a monomial.
             * \return Reference on this.
             */
            Polynomial & operator*=(int numMon);
            
            /**
             * \brief Overload the operator *= to multiply this with a term.
             * \param term: Term.
             * \return Reference on this.
             */
            Polynomial & operator*=(Term<Element> const & term);
            
        private:
            std::forward_list<Term<Element>> _polynomial; /*!< Define a polynomial as a single chained list of Terms. */
    };
    
    /**
     * \brief Overload the operator <<.
     * \return ostream: Stream.
     */
    template <typename Element>
    std::ostream & operator<<(std::ostream & stream, Polynomial<Element> const & polynomial);
}

#include "../src/polynomial.inl"

#endif // F4_POLYNOMIAL_H
