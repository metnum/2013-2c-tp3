/* 
 * File:   MetodoIterativo.h
 * Author: Ramiro
 *
 * Created on 21 de octubre de 2010, 00:58
 */

#ifndef METODOITERATIVO_H
#define	METODOITERATIVO_H

#include <stdexcept>

#include "MatrizRala.h"
#include "CalculadorDeRanking.h"

template <class T>
class MetodoIterativo : public CalculadorDeRanking<T>
{
public:
    /**
     * Constructor.
     * @param unsigned int iteraciones.
     * @param T residuo.
     * @param T uno.
     * @param const MatrizRala<T>* x0 es mejor que este ordenado por columnas.
     */
    MetodoIterativo(unsigned int iteraciones, T residuo, T uno, const MatrizRala<T>* x0)
    {
        this->iteraciones = iteraciones;
        this->residuo = residuo;
        this->uno = uno;
        this->x0 = x0;
    }

    ~MetodoIterativo()
    {
        this->x0 = NULL;
    }

    const MatrizRala<T>* Calcular(const MatrizRala<T>*, const MatrizRala<T>*);

private:
    class Descomposicion
    {
    public:
        const MatrizRala<T> *D_inversa, *L_mas_U;
    };

    Descomposicion* Descomponer(const MatrizRala<T>*);

    unsigned int iteraciones;
    T residuo;
    T uno;
    const MatrizRala<T>* x0;
};

template <class T>
typename MetodoIterativo<T>::Descomposicion* MetodoIterativo<T>::Descomponer(const MatrizRala<T>* A)
{
    typename MatrizRala<T>::Orden* orden_por_filas = new typename MatrizRala<T>::OrdenPorFilas();

    typename MatrizRala<T>::Creador creador_D_inversa(orden_por_filas, A->Filas(), A->Columnas(), A->ElementoNulo());
    typename MatrizRala<T>::Creador creador_L_mas_U(orden_por_filas, A->Filas(), A->Columnas(), A->ElementoNulo());

    typename MatrizRala<T>::Iterador iterador(A);

    //no importa si A esta ordenada por filas o por columnas
    //asi que primero recorro el "orden primario"
    for (iterador.IniciarOrdenPrimario(); iterador.QuedanEnOrdenPrimario(); iterador.SiguienteEnOrdenPrimario())
    {
        //y despues recorro el "orden secundario"
        for (iterador.IniciarOrdenSecundario(); iterador.QuedanEnOrdenSecundario(); iterador.SiguienteEnOrdenSecundario())
        {
            //si es la diagonal invierto el valor y lo guardo en D-1
            if (iterador.FilaActual() == iterador.ColumnaActual())
            {
                creador_D_inversa.Agregar(iterador.FilaActual(), iterador.ColumnaActual(), uno / iterador.ElementoActual());
            }
            //si no guardo el valor en L+U
            else
            {
                creador_L_mas_U.Agregar(iterador.FilaActual(), iterador.ColumnaActual(), iterador.ElementoActual());
            }
        }
    }

    Descomposicion* resultado = new Descomposicion;
    resultado->D_inversa = creador_D_inversa.Crear();
    resultado->L_mas_U = creador_L_mas_U.Crear();

    //cout << "D_inversa" << endl << *(resultado->D_inversa) << endl;
    //cout << "L_mas_U" << endl << *(resultado->L_mas_U) << endl;

    delete orden_por_filas;

    return resultado;
};

/**
 * Calcula el ranking x a partir de A y b resolviendo el sistema Ax = b de modo iterativo.
 * @param const MatrizRala<T>* A es mejor que este ordenada por filas.
 * @param const MatrizRala<T>* b es mejor que este ordenada por filas.
 * @return
 */
template <class T>
const MatrizRala<T>* MetodoIterativo<T>::Calcular(const MatrizRala<T>* A, const MatrizRala<T>* b)
{
    if(A->Columnas() != x0->Filas() || A->Columnas() != b->Filas())
        throw new std::logic_error("Las dimensiones no son adecuadas.");

    if(A->ElementoNulo() != x0->ElementoNulo() || A->ElementoNulo() != b->ElementoNulo())
        throw new std::logic_error("Los elementos nulos no coinciden.");

    //ordenes para las matrices
    typename MatrizRala<T>::Orden* orden_por_filas = new typename MatrizRala<T>::OrdenPorFilas();
    typename MatrizRala<T>::Orden* orden_por_columnas = new typename MatrizRala<T>::OrdenPorColumnas();

    //ordeno b por filas
    const MatrizRala<T> *b_ordenado = b->OrdenadaPorFilas();

    //descompongo A = D + L + U
    //directamente obtengo D inversa y L+U
    Descomposicion* descomposicion = Descomponer(A);

    //variables temporales
    const MatrizRala<T> *x, *tmp1, *tmp2, *xtmp, *Ax, *b_menos_Ax;

    //valor inicial
    x = x0;

    //guarda para para antes de tiempo
    bool continuar = true;

    //residuo en cada iteracion
    T residuo;

    //cantidad de iteraciones ejecutadas
    unsigned int iteracion = 0;
    
    cout.precision(50);

    //por cada iteracion mientras no haya que detenerse
    for (; iteracion < iteraciones && continuar; ++iteracion)
    {
        //tmp1 = (L+U)*x
        //tmp1 es mejor que este ordenada por filas pq se va a restar conta b que esta ordenada por filas
        tmp1 = descomposicion->L_mas_U->MultiplicarPor(x, orden_por_filas);

        //tmp2 = b - tmp1 = b - (L+U)*x
        //tmp2 es mejor que este ordenada por columnas pq se va a multiplicar en la forma A * tmp2
        tmp2 = b_ordenado->Restarle(tmp1, orden_por_columnas);

        //xtmp = D-1 * tmp2 = D-1 * (b - (L+U)*x)
        //tmp2 es mejor que este ordenada por columnas pq se va a convertir en el siguiente x
        xtmp = descomposicion->D_inversa->MultiplicarPor(tmp2, orden_por_columnas);

        if(x != x0) delete x;
        delete tmp1;
        delete tmp2;

        x = xtmp;

        //calculo el residuo
        Ax = A->MultiplicarPor(x, orden_por_filas);
        b_menos_Ax = b_ordenado->Restarle(Ax, orden_por_filas);
        residuo = b_menos_Ax->NormaFrobenius();

        //cout << "iteracion " << i << ": " << " residuo = " << residuo << " / " << this->residuo << endl;

        delete Ax;
        delete b_menos_Ax;

        //verifico si el residuo es suficientemente chico
        continuar = residuo >= this->residuo;

        cout << (iteracion + 1) << " " << residuo << endl;
    }

    //cout << "residuo final = " << residuo << " / " << this->residuo << endl;

    //voy a normalizar (norma 1)
    typename MatrizRala<T>::Creador creador_resultado(orden_por_filas, x->Filas(), 1, x->ElementoNulo());
    const MatrizRala<T>* x_ordenado = x->OrdenadaPorColumnas();
    typename MatrizRala<T>::Iterador iterador(x_ordenado);
    T sumatoria = x->ElementoNulo();

    //busco la suma de todos los elementos
    //itero primero por columnas pq x esta ordenada por columnas
    for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
    {
        for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
        {
            sumatoria+= fabs(iterador.ElementoActual());
        }
    }

    //armo el resultado normalizado
    //itero primero por columnas pq x esta ordenada por columnas
    for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
    {
        for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
        {
            creador_resultado.Agregar(iterador.FilaActual(), iterador.ColumnaActual(), iterador.ElementoActual() / sumatoria);
        }
    }

    const MatrizRala<T>* x_normalizado = creador_resultado.Crear();

    if(x_ordenado != x) delete x_ordenado;
    if(x != x0) delete x;
    if(b_ordenado != b) delete b_ordenado;

    delete orden_por_filas;
    delete orden_por_columnas;
    delete descomposicion;

    return x_normalizado;
}

#endif	/* METODOITERATIVO_H */

