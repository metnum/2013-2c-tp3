/* 
 * File:   MetodoDirecto.h
 * Author: Ramiro
 *
 * Created on 21 de octubre de 2010, 00:52
 */

#ifndef METODODIRECTO_H
#define	METODODIRECTO_H

#include <stdexcept>
#include <vector>

#include "MatrizRala.h"
#include "CalculadorDeRanking.h"

using namespace std;

template <class T>
class MetodoDirecto : public CalculadorDeRanking<T>
{
public:
    const MatrizRala<T>* Calcular(const MatrizRala<T>*, const MatrizRala<T>*);

private:
    class QR
    {
    public:
        QR(const MatrizRala<T>* Q_t, const MatrizRala<T>* R)
        {
            this->Q_t = Q_t;
            this->R = R;
        }

        ~QR()
        {
            delete Q_t;
            delete R;
        }

        const MatrizRala<T>* Q_t;
        const MatrizRala<T>* R;
    };

    const MatrizRala<T>* Givens(unsigned int, unsigned int, const MatrizRala<T>*, typename MatrizRala<T>::Orden*);
    const QR* QRGivens(const MatrizRala<T>*);
    const MatrizRala<T>* Sustitucion(const MatrizRala<T>*, const MatrizRala<T>*, typename MatrizRala<T>::Orden*);
    const MatrizRala<T>* Normalizar(const MatrizRala<T>*, typename MatrizRala<T>::Orden*);
};

template <class T>
const MatrizRala<T>* MetodoDirecto<T>::Sustitucion
(
    const MatrizRala<T>* A,
    const MatrizRala<T>* b,
    typename MatrizRala<T>::Orden* orden
){
    //cantidad de filas
    size_t filas = A->Filas();
    //cantidad de columnas
    size_t columnas = A->Columnas();
    //cantidad de columnas
    T elemento_nulo = A->ElementoNulo();

    //vector incognita
    vector<T>* x = new vector<T>(filas);

    //acumulador temporal
    T acum;

    //recorro desde la ultima fila
    for (int i = filas - 1; i >= 0; i--)
    {
        //si el valor en la diagonal es 0, fijo x como 0 y salteo para no hacer una division por 0
        if (A->ValorEn(i, i) == A->ElementoNulo())
        {
            (*x)[i] = elemento_nulo;
            continue;
        }

        //acum = b[i]
        acum = b->ValorEn(i, 0);

        //para cada elemento j despues de la diagonal
        //acum-= A[i][j] * x[j]
        for (unsigned int j = i + 1; j < columnas; j++)
            acum-= A->ValorEn(i, j) * (*x)[j];

        //finalmente como A[i][i] != 0
        //x[i] = acum / A[i][i]
        (*x)[i] = acum / A->ValorEn(i, i);
    }

    //creador
    typename MatrizRala<T>::Creador creador(orden, filas, 1, elemento_nulo);

    //convierto el vector
    for (size_t i = 0; i < filas; i++)
        creador.Agregar(i, 0, (*x)[i]);

    //libero el vector
    delete x;

    return creador.Crear();
}

template <class T>
const MatrizRala<T>* MetodoDirecto<T>::Normalizar(const MatrizRala<T>* v, typename MatrizRala<T>::Orden* orden)
{
    if(v->Columnas() != 1)
        throw new std::logic_error("La matriz a normalizar debe tener una sola columna.");

    //creador
    typename MatrizRala<T>::Creador creador(orden, v->Filas(), 1, v->ElementoNulo());

    //norma
    T norma = v->ElementoNulo();

    //lo ordeno por filas
    const MatrizRala<T>* ordenado = v->OrdenadaPorFilas();

    //iterador para acceder mas rapido a cada elemento
    typename MatrizRala<T>::Iterador iterador(ordenado);

    //recolecto la norma
    for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
    {
        for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
            norma+= fabs(iterador.ElementoActual());
    }

    //creo el vector normalizado
    for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
    {
        for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
            creador.Agregar(iterador.FilaActual(), 0, iterador.ElementoActual() / norma);
    }

    if(ordenado != v) delete ordenado;

    return creador.Crear();
}

template <class T>
const MatrizRala<T>* MetodoDirecto<T>::Givens
(
    unsigned int i,
    unsigned int j,
    const MatrizRala<T>* A,
    typename MatrizRala<T>::Orden* orden
){
    typename MatrizRala<T>::Creador creador_Q(orden, A->Filas(), A->Columnas(), A->ElementoNulo());

    //Creo identidad
    for(size_t h = 0; h < A->Columnas(); ++h)
    {
        if(h != i && h != j)
            creador_Q.Agregar(h, h, 1);
    }

    T a = A->ValorEn(j, j);
    T b = A->ValorEn(i, j);

    T r = sqrt(a * a + b * b);

    T c = a / r;
    T s = b / r;

    //relleno la matriz
    creador_Q.Agregar(i, i,  c);
    creador_Q.Agregar(j, j,  c);
    creador_Q.Agregar(i, j, -s);
    creador_Q.Agregar(j, i,  s);

    return creador_Q.Crear();
}

template <class T>
const typename MetodoDirecto<T>::QR* MetodoDirecto<T>::QRGivens(const MatrizRala<T>* A)
{
    typename MatrizRala<T>::Orden* orden_por_columnas = new typename MatrizRala<T>::OrdenPorColumnas();
    typename MatrizRala<T>::Orden* orden_por_filas = new typename MatrizRala<T>::OrdenPorFilas();

    //dimensiones
    size_t filas = A->Filas();
    size_t columnas = A->Columnas();
    //elemento nulo
    T elemento_nulo = A->ElementoNulo();

    //matriz R
    //inicialmente es A
    //va ordenada por columnas pq multiplica a derecha
    typename MatrizRala<T>::Creador creador_R_inicial(orden_por_columnas, filas, columnas, elemento_nulo);
    creador_R_inicial.Agregar(*A);
    const MatrizRala<T>* R = creador_R_inicial.Crear();

    //matriz de givens temporal
    const MatrizRala<T>* G;
    //matriz temporal
    const MatrizRala<T>* tmp;

    //matriz Q transpuesta
    //inicialmente es una identidad
    //va ordenada por columnas pq multiplica a derecha
    const MatrizRala<T>* Q_t = identidad(orden_por_columnas, filas, columnas, elemento_nulo);

    //Para cada columna
    for(unsigned int j = 0; j < filas; j++)
    {
    	//Por cada fila
        for(unsigned int i = j + 1; i < columnas; i++)
        {
            //si no vale la pena no hago este paso
            if(R->ValorEn(i, j) == R->ElementoNulo())
                continue;

            //obtengo la matriz de givens
            //va ordenada por filas pq multiplica a izquierda
            G = Givens(i, j, R, orden_por_filas);

            //obtengo la nueva R multiplicando G * R
            //va ordenada por columnas pq multiplica a derecha
            tmp = G->MultiplicarPor(R, orden_por_columnas);
            delete R;
            R = tmp;

            //copio R pero le saco la posicion i j
            typename MatrizRala<T>::Creador creador_R(orden_por_columnas, filas, columnas, elemento_nulo);
            creador_R.Agregar(*R);
            creador_R.Sacar(i, j);
            tmp = creador_R.Crear();
            delete R;
            R = tmp;

            //acumulo las givens para formar Q_t
            //va ordenada por columnas pq multiplica a derecha
            tmp = G->MultiplicarPor(Q_t, orden_por_columnas);
            delete Q_t;
            Q_t = tmp;

            //libero memoria
            delete G;
        }
    }

    //libero memoria
    delete orden_por_columnas;

    return new QR(Q_t, R);
}

template <class T>
const MatrizRala<T>* MetodoDirecto<T>::Calcular(const MatrizRala<T>* A, const MatrizRala<T>* b)
{
    if(A->Columnas() != b->Filas())
        throw new std::logic_error("Las dimensiones no son adecuadas.");

    if(A->ElementoNulo() != b->ElementoNulo())
        throw new std::logic_error("Los elementos nulos no coinciden.");

    typename MatrizRala<T>::Orden* orden_por_columnas = new typename MatrizRala<T>::OrdenPorColumnas();
    typename MatrizRala<T>::Orden* orden_por_filas = new typename MatrizRala<T>::OrdenPorFilas();

    /*
     * calcula Q_t y R tal que:
     *
     * 1) A = Q * R
     * 2) R es triangular superior
     * 3) Q es ortogonal
     * 
     * Ademas si A * x = b entonces:
     * 
     * Q_t * A * x = Q_t * b
     * 
     * luego como A = Q * R y Q_t * Q = I entonces:
     * 
     * R * x = Q_t * b
     */
    const QR* QR = QRGivens(A);

    //cout << "Q_t" << endl << *(QR->Q_t) << endl;

    //cout << "R" << endl << *(QR->R) << endl;

    //calculo Q_t * b
    const MatrizRala<T>* Q_t_b = QR->Q_t->MultiplicarPor(b, orden_por_filas);

    //resuelvo el sistema R * x = Q_t * b
    //simplemente hago sustitucion hacia atras
    const MatrizRala<T>* x = Sustitucion(QR->R, Q_t_b, orden_por_filas);

    //cout << "x (sin normalizar)" << endl << *x << endl;

    //normalizo x
    const MatrizRala<T>* x_normalizado = Normalizar(x, orden_por_columnas);

    delete orden_por_columnas;
    delete orden_por_filas;
    delete QR;
    delete Q_t_b;
    delete x;

    return x_normalizado;
}

#endif	/* METODODIRECTO_H */
