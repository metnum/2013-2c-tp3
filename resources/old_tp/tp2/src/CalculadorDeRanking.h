/* 
 * File:   CalculadorDeRanking.h
 * Author: Ramiro
 *
 * Created on 21 de octubre de 2010, 00:50
 */

#ifndef CALCULADORDERANKING_H
#define	CALCULADORDERANKING_H

#include "MatrizRala.h"

using namespace std;

template <class T>
class CalculadorDeRanking
{
public:
    virtual const MatrizRala<T>* Calcular(const MatrizRala<T>*, const MatrizRala<T>*) = 0;
};

#endif	/* CALCULADORDERANKING_H */

