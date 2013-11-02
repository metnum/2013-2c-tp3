/* 
 * File:   MatrizRala.h
 * Author: ramiro
 *
 * Created on 4 de octubre de 2010, 16:17
 */

#ifndef _MATRIZRALA_H
#define	_MATRIZRALA_H

#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <stdexcept>

using namespace std;

template <class T> class MatrizRala
{
public:
    class Orden;
    class Posicion;
    class Creador;
    class Iterador;

    MatrizRala
    (
        Orden* orden,
        T* elementos,
        unsigned int* inicios_orden_primario,
        unsigned int* indices_orden_primario,
        unsigned int* indices_orden_secundario,
        T elemento_nulo,
        size_t filas,
        size_t columnas,
        size_t longitud_orden_primario,
        size_t longitud
    ){
        this->orden = orden;
        this->elementos = elementos;
        this->inicios_orden_primario = inicios_orden_primario;
        this->indices_orden_primario = indices_orden_primario;
        this->indices_orden_secundario = indices_orden_secundario;
        this->elemento_nulo = elemento_nulo;
        this->filas = filas;
        this->columnas = columnas;
        this->longitud_orden_primario = longitud_orden_primario;
        this->longitud = longitud;
    }

    ~MatrizRala()
    {
        delete this->orden;
        delete[] this->elementos;
        delete[] this->inicios_orden_primario;
        delete[] this->indices_orden_primario;
        delete[] this->indices_orden_secundario;

        this->orden = NULL;
        this->elementos = NULL;
        this->inicios_orden_primario = NULL;
        this->indices_orden_primario = NULL;
        this->indices_orden_secundario = NULL;
    }

    const MatrizRala<T>* SumarA(const MatrizRala<T>*, Orden*) const;
    const MatrizRala<T>* SumarA(const MatrizRala<T>&, Orden*) const;
    const MatrizRala<T>* Restarle(const MatrizRala<T>*, Orden*) const;
    const MatrizRala<T>* Restarle(const MatrizRala<T>&, Orden*) const;

    const MatrizRala<T>* MultiplicarPor(const MatrizRala<T>*, Orden*) const;
    const MatrizRala<T>* MultiplicarPor(const MatrizRala<T>&, Orden*) const;
    const MatrizRala<T>* MultiplicarPor(const T&, Orden*) const;
    
    vector<T>* MultiplicarPor(const vector<T>&) const;
    
    template <class T2>
    friend vector<T2>* Multiplicar(const vector<T2>&, const MatrizRala<T2>&);
    template <class T2>
    friend vector<T2>* Multiplicar(const vector<T2>&, const MatrizRala<T2>*);

    template <class T2>
    friend const MatrizRala<T2>* Multiplicar(const T2&, const MatrizRala<T2>&, typename MatrizRala<T2>::Orden*);
    template <class T2>
    friend const MatrizRala<T2>* Multiplicar(const T2&, const MatrizRala<T2>*, typename MatrizRala<T2>::Orden*);

    template <class T2>
    friend ostream& operator<<(ostream&, const MatrizRala<T2>&);

    T ValorEn(size_t, size_t) const;
    T operator()(size_t, size_t) const;

    const MatrizRala<T>* OrdenadaPorFilas() const;
    const MatrizRala<T>* OrdenadaPorColumnas() const;
    const MatrizRala<T>* Ordenar(Orden*) const;

    const MatrizRala<T>* Trasponer(Orden*) const;

    T NormaFrobenius() const;

    size_t Longitud() const { return longitud; }

    size_t Filas() const { return filas; }
    size_t Columnas() const { return columnas; }
    T ElementoNulo() const { return elemento_nulo; }

    class Orden
    {
    public:
        virtual Orden* duplicar() = 0;

        virtual bool PrimeraEsMenor(const Posicion&, const Posicion&) = 0;
        
        virtual size_t OrdenPrimario(const Posicion&) = 0;
        virtual size_t OrdenSecundario(const Posicion&) = 0;

        virtual void IniciarFilas(Iterador*) = 0;
        virtual void IniciarColumnas(Iterador*) = 0;
        
        virtual void SiguienteFila(Iterador*) = 0;
        virtual void SiguienteColumna(Iterador*) = 0;

        virtual bool QuedanFilas(Iterador*) = 0;
        virtual bool QuedanColumnas(Iterador*) = 0;

        virtual unsigned int FilaActual(Iterador*) = 0;
        virtual unsigned int ColumnaActual(Iterador*) = 0;

        virtual const MatrizRala<T>* OrdenarPorFilas(const MatrizRala<T>*) = 0;
        virtual const MatrizRala<T>* OrdenarPorColumnas(const MatrizRala<T>*) = 0;
    };

    class OrdenPorFilas : public Orden
    {
    public:
        Orden* duplicar() { return new OrdenPorFilas(); }

        bool PrimeraEsMenor(const Posicion&, const Posicion&);

        size_t OrdenPrimario(const Posicion&);
        size_t OrdenSecundario(const Posicion&);

        void IniciarFilas(Iterador*);
        void IniciarColumnas(Iterador*);

        void SiguienteFila(Iterador*);
        void SiguienteColumna(Iterador*);

        bool QuedanFilas(Iterador*);
        bool QuedanColumnas(Iterador*);

        unsigned int FilaActual(Iterador*);
        unsigned int ColumnaActual(Iterador*);

        const MatrizRala<T>* OrdenarPorFilas(const MatrizRala<T>*);
        const MatrizRala<T>* OrdenarPorColumnas(const MatrizRala<T>*);
    };

    class OrdenPorColumnas : public Orden
    {
    public:
        Orden* duplicar() { return new OrdenPorColumnas(); }

        bool PrimeraEsMenor(const Posicion&, const Posicion&);

        size_t OrdenPrimario(const Posicion&);
        size_t OrdenSecundario(const Posicion&);

        void IniciarFilas(Iterador*);
        void IniciarColumnas(Iterador*);

        void SiguienteFila(Iterador*);
        void SiguienteColumna(Iterador*);

        bool QuedanFilas(Iterador*);
        bool QuedanColumnas(Iterador*);

        unsigned int FilaActual(Iterador*);
        unsigned int ColumnaActual(Iterador*);

        const MatrizRala<T>* OrdenarPorFilas(const MatrizRala<T>*);
        const MatrizRala<T>* OrdenarPorColumnas(const MatrizRala<T>*);
    };

    class Posicion
    {
    public:
        Posicion(size_t fila, size_t columna, Orden* orden)
        {
            this->fila = fila;
            this->columna = columna;
            this->orden = orden;
        }

        bool operator<(const Posicion&) const;

        size_t Fila() const { return fila; }
        size_t Columna() const { return columna; }

    private:
        Orden* orden;

        size_t fila;
        size_t columna;
    };

    class Creador
    {
    public:
        typedef map<Posicion, T> Diccionario;
        typedef typename Diccionario::iterator IteradorInterno;

        Creador(Orden* orden, size_t filas, size_t columnas, T elemento_nulo)
        {
            this->orden = orden;
            this->filas = filas;
            this->columnas = columnas;
            this->elemento_nulo = elemento_nulo;
        }

        ~Creador()
        {
            this->orden = NULL;
        }

        void Agregar(size_t fila, size_t columna, T elemento);
        void Agregar(const MatrizRala<T>&);
        void Sacar(size_t fila, size_t columna);

        const MatrizRala<T>* Crear();
    private:
        Orden* orden;
        Diccionario diccionario;
        size_t filas;
        size_t columnas;
        T elemento_nulo;
    };

    class Iterador
    {
    public:
        Iterador(const MatrizRala<T>* matriz)
        {
            this->matriz = matriz;
        }

        void IniciarFilas();
        void IniciarColumnas();
        
        void SiguienteFila();
        void SiguienteColumna();

        bool QuedanFilas();
        bool QuedanColumnas();

        unsigned int FilaActual();
        unsigned int ColumnaActual();

        T ElementoActual();

        void IniciarOrdenPrimario();
        void IniciarOrdenSecundario();
        void SiguienteEnOrdenPrimario();
        void SiguienteEnOrdenSecundario();

        bool QuedanEnOrdenPrimario();
        bool QuedanEnOrdenSecundario();

        unsigned int OrdenPrimarioActual();
        unsigned int OrdenSecundarioActual();

    private:
        const MatrizRala<T>* matriz;

        unsigned int indice_orden_primario;
        unsigned int inicio_orden_primario;
        unsigned int fin_orden_primario;
        unsigned int indice_orden_secundario;
    };

private:
    class OperacionBinaria
    {
    public:
        virtual T operator()(const T&, const T&) const = 0;
    };

    class Suma : public OperacionBinaria
    {
    public:
        T operator()(const T& t1, const T& t2) const { return t1 + t2; }
    };

    class Resta : public OperacionBinaria
    {
    public:
        T operator()(const T& t1, const T& t2) const { return t1 - t2; }
    };

    const MatrizRala<T>* OperacionBinariaCon(const OperacionBinaria&, const MatrizRala<T>&, Orden*) const;

    Orden* orden;
    T* elementos;
    unsigned int* inicios_orden_primario;
    unsigned int* indices_orden_primario;
    unsigned int* indices_orden_secundario;
    T elemento_nulo;
    size_t filas;
    size_t columnas;
    size_t longitud_orden_primario;
    size_t longitud;
};

template <class T>
const MatrizRala<T>* identidad(typename MatrizRala<T>::Orden* orden, size_t filas, size_t columnas, T elemento_nulo)
{
    typename MatrizRala<T>::Creador creador(orden, filas, columnas, elemento_nulo);

    for(size_t i = 0; i < filas; ++i)
        creador.Agregar(i, i, 1.0);

    return creador.Crear();
}

template <class T>
ostream& operator<<(ostream& os, const MatrizRala<T>& matriz)
{
    size_t filas = matriz.Filas();
    size_t columnas = matriz.Columnas();

    for(unsigned int i = 0; i < filas; i++)
    {
        if(i > 0)
            os << endl;

        if(columnas > 0)
            os << matriz.ValorEn(i, 0);

        for(unsigned int j = 1; j < columnas; j++)
        {
            os << " " << matriz.ValorEn(i, j);
        }
    }

    return os;
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::SumarA(const MatrizRala<T>* otra, Orden* orden_resultado) const
{
    return SumarA(*otra, orden_resultado);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::SumarA(const MatrizRala<T>& otra, Orden* orden_resultado) const
{
    return OperacionBinariaCon(Suma(), otra, orden_resultado);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::Restarle(const MatrizRala<T>* otra, Orden* orden_resultado) const
{
    return Restarle(*otra, orden_resultado);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::Restarle(const MatrizRala<T>& otra, Orden* orden_resultado) const
{
    return OperacionBinariaCon(Resta(), otra, orden_resultado);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OperacionBinariaCon(const OperacionBinaria& operacion_binaria, const MatrizRala<T>& otra, Orden* orden_resultado) const
{
    if(filas != otra.filas)
        throw new length_error("Cantidad de filas incompatibles.");

    if(columnas != otra.columnas)
        throw new length_error("Cantidad de columnas incompatibles.");

    if(elemento_nulo != otra.elemento_nulo)
        throw new logic_error("Elementos nulos incompatibles.");

    Creador creador(orden_resultado, filas, columnas, elemento_nulo);

    //ordeno ambas por filas
    const MatrizRala<T>* A = this->OrdenadaPorFilas();
    const MatrizRala<T>* B = otra.OrdenadaPorFilas();

    //iteradores
    Iterador iterador_A(A);
    Iterador iterador_B(B);

    //comienzo la iteracion por filas
    iterador_A.IniciarFilas();
    iterador_B.IniciarFilas();

    //mientras que ambas tengan mas filas
    while(iterador_A.QuedanFilas() && iterador_B.QuedanFilas())
    {
        //si estan en la misma fila
        if(iterador_A.FilaActual() == iterador_B.FilaActual())
        {
            //inicio las columnas
            iterador_A.IniciarColumnas();
            iterador_B.IniciarColumnas();

            //mientras que ambas tengan mas columnas
            while(iterador_A.QuedanColumnas() && iterador_B.QuedanColumnas())
            {
                //si estan en la misma columna
                if(iterador_A.ColumnaActual() == iterador_B.ColumnaActual())
                {
                    //resultado[i][j] = operacion_binaria(A[i][j], B[i][j])
                    creador.Agregar
                    (
                        iterador_A.FilaActual(),
                        iterador_A.ColumnaActual(),
                        operacion_binaria(iterador_A.ElementoActual(), iterador_B.ElementoActual())
                    );

                    //avanzo las columnas
                    iterador_A.SiguienteColumna();
                    iterador_B.SiguienteColumna();
                }
                //si A esta en una columna anterior
                else if(iterador_A.ColumnaActual() < iterador_B.ColumnaActual())
                {
                    //resultado[i][j] = operacion_binaria(A[i][j], 0)
                    creador.Agregar
                    (
                        iterador_A.FilaActual(),
                        iterador_A.ColumnaActual(),
                        operacion_binaria(iterador_A.ElementoActual(), elemento_nulo)
                    );

                    //avanzo la columna de A
                    iterador_A.SiguienteColumna();
                }
                //si B esta en una columna anterior
                else if(iterador_B.ColumnaActual() < iterador_A.ColumnaActual())
                {
                    //resultado[i][j] = operacion_binaria(0, B[i][j])
                    creador.Agregar
                    (
                        iterador_B.FilaActual(),
                        iterador_B.ColumnaActual(),
                        operacion_binaria(elemento_nulo, iterador_B.ElementoActual())
                    );

                    //avanzo la columna de B
                    iterador_B.SiguienteColumna();
                }
            }

            //agrego si sobran columnas de A
            for(; iterador_A.QuedanColumnas(); iterador_A.SiguienteColumna())
            {
                //resultado[i][j] = operacion_binaria(A[i][j], 0)
                creador.Agregar
                (
                    iterador_A.FilaActual(),
                    iterador_A.ColumnaActual(),
                    operacion_binaria(iterador_A.ElementoActual(), elemento_nulo)
                );
            }

            //agrego si sobran columnas de B
            for(; iterador_B.QuedanColumnas(); iterador_B.SiguienteColumna())
            {
                //resultado[i][j] = operacion_binaria(0, B[i][j])
                creador.Agregar
                (
                    iterador_B.FilaActual(),
                    iterador_B.ColumnaActual(),
                    operacion_binaria(elemento_nulo, iterador_B.ElementoActual())
                );
            }

            //avanzo las filas
            iterador_A.SiguienteFila();
            iterador_B.SiguienteFila();
        }
        //si A esta en una fila anterior
        else if(iterador_A.FilaActual() < iterador_B.FilaActual())
        {
            for(iterador_A.IniciarColumnas(); iterador_A.QuedanColumnas(); iterador_A.SiguienteColumna())
            {
                //resultado[i][j] = operacion_binaria(A[i][j], 0)
                creador.Agregar
                (
                    iterador_A.FilaActual(),
                    iterador_A.ColumnaActual(),
                    operacion_binaria(iterador_A.ElementoActual(), elemento_nulo)
                );
            }

            //avanzo la fila de A
            iterador_A.SiguienteFila();
        }
        //si B esta en una fila anterior
        else if(iterador_B.FilaActual() < iterador_A.FilaActual())
        {
            for(iterador_B.IniciarColumnas(); iterador_B.QuedanColumnas(); iterador_B.SiguienteColumna())
            {
                //resultado[i][j] = operacion_binaria(0, B[i][j])
                creador.Agregar
                (
                    iterador_B.FilaActual(),
                    iterador_B.ColumnaActual(),
                    operacion_binaria(elemento_nulo, iterador_B.ElementoActual())
                );
            }

            //avanzo la fila de B
            iterador_B.SiguienteFila();
        }
    }

    //agrego si sobran filas de A
    for(; iterador_A.QuedanFilas(); iterador_A.SiguienteFila())
    {
        for(iterador_A.IniciarColumnas(); iterador_A.QuedanColumnas(); iterador_A.SiguienteColumna())
        {
            //resultado[i][j] = operacion_binaria(A[i][j], 0)
            creador.Agregar
            (
                iterador_A.FilaActual(),
                iterador_A.ColumnaActual(),
                operacion_binaria(iterador_A.ElementoActual(), elemento_nulo)
            );
        }
    }

    //agrego si sobran filas de B
    for(; iterador_B.QuedanFilas(); iterador_B.SiguienteFila())
    {
        for(iterador_B.IniciarColumnas(); iterador_B.QuedanColumnas(); iterador_B.SiguienteColumna())
        {
            //resultado[i][j] = operacion_binaria(0, B[i][j])
            creador.Agregar
            (
                iterador_B.FilaActual(),
                iterador_B.ColumnaActual(),
                operacion_binaria(elemento_nulo, iterador_B.ElementoActual())
            );
        }
    }

    if(A != this) delete A;
    if(B != &otra) delete B;

    return creador.Crear();
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::MultiplicarPor(const MatrizRala<T>* otra, Orden* orden_resultado) const
{
    return MultiplicarPor(*otra, orden_resultado);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::MultiplicarPor(const MatrizRala<T>& otra, Orden* orden_resultado) const
{
    if(columnas != otra.filas)
        throw new length_error("Cantidad de filas y columnas incompatibles.");

    if(elemento_nulo != otra.elemento_nulo)
        throw new logic_error("Elementos nulos incompatibles.");

    const MatrizRala<T>* A = this->OrdenadaPorFilas();
    const MatrizRala<T>* B = otra.OrdenadaPorColumnas();
    
    Creador creador(orden_resultado, filas, otra.columnas, elemento_nulo);

    Iterador iterador_a(A);
    Iterador iterador_b(B);

    T valor;

    for(iterador_a.IniciarFilas(); iterador_a.QuedanFilas(); iterador_a.SiguienteFila())
    {
        for(iterador_b.IniciarColumnas(); iterador_b.QuedanColumnas(); iterador_b.SiguienteColumna())
        {
            iterador_a.IniciarColumnas();
            iterador_b.IniciarFilas();

            valor = 0;

            while(iterador_a.QuedanColumnas() && iterador_b.QuedanFilas())
            {
                if(iterador_a.ColumnaActual() == iterador_b.FilaActual())
                {
                    valor+= iterador_a.ElementoActual() * iterador_b.ElementoActual();

                    iterador_a.SiguienteColumna();
                    iterador_b.SiguienteFila();
                }
                else if(iterador_a.ColumnaActual() < iterador_b.FilaActual())
                {
                    iterador_a.SiguienteColumna();
                }
                else
                {
                    iterador_b.SiguienteFila();
                }
            }

            creador.Agregar(iterador_a.FilaActual(), iterador_b.ColumnaActual(), valor);
        }
    }

    if(A != this)
        delete A;

    if(B != &otra)
        delete B;

    return creador.Crear();
}

template <class T>
vector<T>* MatrizRala<T>::MultiplicarPor(const vector<T>& x) const
{
    if(columnas != x.size())
        throw new length_error("Cantidad de filas y columnas incompatibles.");

    vector<T>* b = new vector<T>(filas);

    for(unsigned int fila = 0; fila < filas; fila++)
        (*b)[fila] = elemento_nulo;

    const MatrizRala<T>* ordenada = this->OrdenadaPorFilas();
    Iterador iterador(ordenada);
    T valor;

    for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
    {
        valor = elemento_nulo;

        for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
        {
            valor+= iterador.ElementoActual() * x[iterador.ColumnaActual()];
        }

        (*b)[iterador.FilaActual()] = valor;
    }

    if(ordenada != this)
        delete ordenada;

    return b;
}

template <class T>
vector<T>* Multiplicar(const vector<T>& x, const MatrizRala<T>* A)
{
    return Multiplicar(x, *A);
}

template <class T>
vector<T>* Multiplicar(const vector<T>& x, const MatrizRala<T>& A)
{
    vector<T>* b = new vector<T>(A.columnas);

    const MatrizRala<T>* ordenada = A.OrdenadaPorColumnas();
    typename MatrizRala<T>::Iterador iterador(ordenada);
    T valor;

    for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
    {
        valor = A.elemento_nulo;

        for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
        {
            valor+= iterador.ElementoActual() * x[iterador.FilaActual()];
        }

        (*b)[iterador.ColumnaActual()] = valor;
    }

    if(ordenada != &A)
        delete ordenada;

    return b;
}

template <class T>
const MatrizRala<T>* Multiplicar(const T& c, const MatrizRala<T>* A, typename MatrizRala<T>::Orden* orden_resultado)
{
    return Multiplicar(c, *A, orden_resultado);
}

template <class T>
const MatrizRala<T>* Multiplicar(const T& c, const MatrizRala<T>& A, typename MatrizRala<T>::Orden* orden_resultado)
{
    return A.MultiplicarPor(c, orden_resultado);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::MultiplicarPor(const T& c, Orden* orden_resultado) const
{
    Creador creador(orden_resultado, filas, columnas, elemento_nulo);

    Iterador iterador(this);

    for(iterador.IniciarOrdenPrimario(); iterador.QuedanEnOrdenPrimario(); iterador.SiguienteEnOrdenPrimario())
    {
        for(iterador.IniciarOrdenSecundario(); iterador.QuedanEnOrdenSecundario(); iterador.SiguienteEnOrdenSecundario())
        {
            creador.Agregar(iterador.FilaActual(), iterador.ColumnaActual(), iterador.ElementoActual() * c);
        }
    }

    return creador.Crear();
}

template <class T>
T MatrizRala<T>::ValorEn(size_t fila, size_t columna) const
{
    if(longitud == 0) return elemento_nulo;

    Posicion p(fila, columna, orden);
    T primario = orden->OrdenPrimario(p);
    T secundario = orden->OrdenSecundario(p);

    bool encontrado = false;

    size_t inicio_primario = 0;
    size_t fin_primario = longitud_orden_primario - 1;
    size_t mitad_primario;
    size_t indice_primario = 0;

    while
    (
        inicio_primario <= fin_primario
        && !encontrado
        && indices_orden_primario[inicio_primario] <= primario
        && primario <= indices_orden_primario[fin_primario]
    ){
        if(indices_orden_primario[inicio_primario] == primario)
        {
            indice_primario = inicio_primario;
            encontrado = true;
        }
        else if(indices_orden_primario[fin_primario] == primario)
        {
            indice_primario = fin_primario;
            encontrado = true;
        }
        else
        {
            mitad_primario = floor((inicio_primario + fin_primario) / 2.0f);

            if(indices_orden_primario[mitad_primario] == primario)
            {
                indice_primario = mitad_primario;
                encontrado = true;
            }
            else if(indices_orden_primario[mitad_primario] < primario)
                inicio_primario = mitad_primario + 1;
            else
                fin_primario = mitad_primario;
        }
    }

    if(encontrado)
    {
        encontrado = false;

        size_t inicio_secundario = inicios_orden_primario[indice_primario];
        size_t fin_secundario = inicios_orden_primario[indice_primario + 1] - 1;
        size_t mitad_secundario;
        size_t indice_secundario = 0;

        while
        (
            inicio_secundario <= fin_secundario
            && !encontrado
            && indices_orden_secundario[inicio_secundario] <= secundario
            && secundario <= indices_orden_secundario[fin_secundario]
        ){
            if(indices_orden_secundario[inicio_secundario] == secundario)
            {
                indice_secundario = inicio_secundario;
                encontrado = true;
            }
            else if(indices_orden_secundario[fin_secundario] == secundario)
            {
                indice_secundario = fin_secundario;
                encontrado = true;
            }
            else
            {
                mitad_secundario = floor((inicio_secundario + fin_secundario) / 2.0f);

                if(indices_orden_secundario[mitad_secundario] == secundario)
                {
                    indice_secundario = mitad_secundario;
                    encontrado = true;
                }
                else if(indices_orden_secundario[mitad_secundario] < secundario)
                    inicio_secundario = mitad_secundario + 1;
                else
                    fin_secundario = mitad_secundario;
            }
        }

        if(encontrado) return elementos[indice_secundario];
    }

    return elemento_nulo;
}

template <class T>
T MatrizRala<T>::operator ()(size_t fila, size_t columna) const
{
    return ValorEn(fila, columna);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OrdenadaPorFilas() const
{
    return orden->OrdenarPorFilas(this);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OrdenadaPorColumnas() const
{
    return orden->OrdenarPorColumnas(this);
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::Ordenar(Orden* orden) const
{
    Creador creador(orden, filas, columnas, elemento_nulo);
    creador.Agregar(*this);
    return creador.Crear();
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::Trasponer(Orden* orden) const
{
    Creador creador(orden, columnas, filas, elemento_nulo);

    Iterador iterador(this);

    for(iterador.IniciarOrdenPrimario(); iterador.QuedanEnOrdenPrimario(); iterador.SiguienteEnOrdenPrimario())
    {
        for(iterador.IniciarOrdenSecundario(); iterador.QuedanEnOrdenSecundario(); iterador.SiguienteEnOrdenSecundario())
        {
            creador.Agregar(iterador.ColumnaActual(), iterador.FilaActual(), iterador.ElementoActual());
        }
    }

    return creador.Crear();
}

template <class T>
T MatrizRala<T>::NormaFrobenius() const
{
    T sumatoria = elemento_nulo;
    T elemento_actual;

    Iterador iterador(this);

    for(iterador.IniciarOrdenPrimario(); iterador.QuedanEnOrdenPrimario(); iterador.SiguienteEnOrdenPrimario())
    {
        for(iterador.IniciarOrdenSecundario(); iterador.QuedanEnOrdenSecundario(); iterador.SiguienteEnOrdenSecundario())
        {
            elemento_actual = iterador.ElementoActual();
            sumatoria+= elemento_actual * elemento_actual;
        }
    }

    return sqrt(sumatoria);
}

template <class T>
bool MatrizRala<T>::OrdenPorFilas::PrimeraEsMenor(const Posicion& p1, const Posicion& p2)
{
    return p1.Fila() < p2.Fila() || (p1.Fila() == p2.Fila() && p1.Columna() < p2.Columna());
}

template <class T>
size_t MatrizRala<T>::OrdenPorFilas::OrdenPrimario(const Posicion& p)
{
    return p.Fila();
}

template <class T>
size_t MatrizRala<T>::OrdenPorFilas::OrdenSecundario(const Posicion& p)
{
    return p.Columna();
}

template <class T>
void MatrizRala<T>::OrdenPorFilas::IniciarFilas(Iterador* iterador)
{
    iterador->IniciarOrdenPrimario();
}

template <class T>
void MatrizRala<T>::OrdenPorFilas::IniciarColumnas(Iterador* iterador)
{
    iterador->IniciarOrdenSecundario();
}

template <class T>
void MatrizRala<T>::OrdenPorFilas::SiguienteFila(Iterador* iterador)
{
    iterador->SiguienteEnOrdenPrimario();
}

template <class T>
void MatrizRala<T>::OrdenPorFilas::SiguienteColumna(Iterador* iterador)
{
    iterador->SiguienteEnOrdenSecundario();
}

template <class T>
bool MatrizRala<T>::OrdenPorFilas::QuedanFilas(Iterador* iterador)
{
    return iterador->QuedanEnOrdenPrimario();
}

template <class T>
bool MatrizRala<T>::OrdenPorFilas::QuedanColumnas(Iterador* iterador)
{
    return iterador->QuedanEnOrdenSecundario();
}

template <class T>
unsigned int MatrizRala<T>::OrdenPorFilas::FilaActual(Iterador* iterador)
{
    return iterador->OrdenPrimarioActual();
}

template <class T>
unsigned int MatrizRala<T>::OrdenPorFilas::ColumnaActual(Iterador* iterador)
{
    return iterador->OrdenSecundarioActual();
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OrdenPorFilas::OrdenarPorFilas(const MatrizRala<T>* matriz)
{
    return matriz;
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OrdenPorFilas::OrdenarPorColumnas(const MatrizRala<T>* matriz)
{
    cerr << "Atencion: cambiando orden de filas a columnas." << endl;

    Orden* orden_opuesto = new OrdenPorColumnas();
    const MatrizRala<T>* matriz_ordenada_al_revez = matriz->Ordenar(orden_opuesto);
    delete orden_opuesto;
    return matriz_ordenada_al_revez;
}

template <class T>
bool MatrizRala<T>::OrdenPorColumnas::PrimeraEsMenor(const Posicion& p1, const Posicion& p2)
{
    return p1.Columna() < p2.Columna() || (p1.Columna() == p2.Columna() && p1.Fila() < p2.Fila());
}

template <class T>
size_t MatrizRala<T>::OrdenPorColumnas::OrdenPrimario(const Posicion& p)
{
    return p.Columna();
}

template <class T>
size_t MatrizRala<T>::OrdenPorColumnas::OrdenSecundario(const Posicion& p)
{
    return p.Fila();
}

template <class T>
void MatrizRala<T>::OrdenPorColumnas::IniciarFilas(Iterador* iterador)
{
    iterador->IniciarOrdenSecundario();
}

template <class T>
void MatrizRala<T>::OrdenPorColumnas::IniciarColumnas(Iterador* iterador)
{
    iterador->IniciarOrdenPrimario();
}

template <class T>
void MatrizRala<T>::OrdenPorColumnas::SiguienteFila(Iterador* iterador)
{
    iterador->SiguienteEnOrdenSecundario();
}

template <class T>
void MatrizRala<T>::OrdenPorColumnas::SiguienteColumna(Iterador* iterador)
{
    iterador->SiguienteEnOrdenPrimario();
}

template <class T>
bool MatrizRala<T>::OrdenPorColumnas::QuedanFilas(Iterador* iterador)
{
    return iterador->QuedanEnOrdenSecundario();
}

template <class T>
bool MatrizRala<T>::OrdenPorColumnas::QuedanColumnas(Iterador* iterador)
{
    return iterador->QuedanEnOrdenPrimario();
}

template <class T>
unsigned int MatrizRala<T>::OrdenPorColumnas::FilaActual(Iterador* iterador)
{
    return iterador->OrdenSecundarioActual();
}

template <class T>
unsigned int MatrizRala<T>::OrdenPorColumnas::ColumnaActual(Iterador* iterador)
{
    return iterador->OrdenPrimarioActual();
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OrdenPorColumnas::OrdenarPorFilas(const MatrizRala<T>* matriz)
{
    cerr << "Atencion: cambiando orden de columnas a filas." << endl;

    Orden* orden_opuesto = new OrdenPorFilas();
    const MatrizRala<T>* matriz_ordenada_al_revez = matriz->Ordenar(orden_opuesto);
    delete orden_opuesto;
    return matriz_ordenada_al_revez;
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::OrdenPorColumnas::OrdenarPorColumnas(const MatrizRala<T>* matriz)
{
    return matriz;
}

template <class T>
bool MatrizRala<T>::Posicion::operator<(const Posicion& otra) const
{
    return orden->PrimeraEsMenor(*this, otra);
}

template <class T>
void MatrizRala<T>::Creador::Agregar(size_t fila, size_t columna, T elemento)
{
    if(fila >= filas)
        throw new out_of_range("Fila fuera de rango: " + fila);

    if(columna >= columnas)
        throw new out_of_range("Columna fuera de rango: " + columna);

    if(elemento != elemento_nulo)
        diccionario[Posicion(fila, columna, orden)] = elemento;
}

template <class T>
void MatrizRala<T>::Creador::Agregar(const MatrizRala<T>& otra)
{
    if(otra.filas > filas)
        throw new out_of_range("La matriz que se desea agregar tiene mas filas de lo que se esperaba.");

    if(otra.columnas > columnas)
        throw new out_of_range("La matriz que se desea agregar tiene mas columnas de lo que se esperaba.");

    Iterador iterador(&otra);

    for(iterador.IniciarOrdenPrimario(); iterador.QuedanEnOrdenPrimario(); iterador.SiguienteEnOrdenPrimario())
    {
        for(iterador.IniciarOrdenSecundario(); iterador.QuedanEnOrdenSecundario(); iterador.SiguienteEnOrdenSecundario())
        {
            Agregar(iterador.FilaActual(), iterador.ColumnaActual(), iterador.ElementoActual());
        }
    }
}

template <class T>
void MatrizRala<T>::Creador::Sacar(size_t fila, size_t columna)
{
    diccionario.erase(Posicion(fila, columna, orden));
}

template <class T>
const MatrizRala<T>* MatrizRala<T>::Creador::Crear()
{
    T* elementos = new T[diccionario.size()];
    unsigned int* indices_orden_secundario = new unsigned int[diccionario.size()];

    size_t dimension_maxima = max(filas, columnas);

    unsigned int* inicios_orden_primario_temporal = new unsigned int[dimension_maxima];
    unsigned int* indices_orden_primario_temporal = new unsigned int[dimension_maxima];

    IteradorInterno iterador;

    unsigned int indice_elemento = 0;

    size_t actual;
    size_t longitud = 0;
    size_t anterior = 0;

    bool primera = true;

    for(iterador = diccionario.begin(); iterador != diccionario.end(); iterador++)
    {
        actual = orden->OrdenPrimario((*iterador).first);

        elementos[indice_elemento] = (*iterador).second;
        indices_orden_secundario[indice_elemento] = orden->OrdenSecundario((*iterador).first);

        if(primera || anterior != actual)
        {
            primera = false;

            anterior = actual;

            inicios_orden_primario_temporal[longitud] = indice_elemento;
            indices_orden_primario_temporal[longitud] = actual;
            
            longitud++;
        }

        indice_elemento++;
    }

    unsigned int* inicios_orden_primario = new unsigned int[longitud + 1];
    unsigned int* indices_orden_primario = new unsigned int[longitud];

    for(actual = 0; actual < longitud; actual++)
    {
        inicios_orden_primario[actual] = inicios_orden_primario_temporal[actual];
        indices_orden_primario[actual] = indices_orden_primario_temporal[actual];
    }

    inicios_orden_primario[longitud] = diccionario.size();

    return new MatrizRala<T>(
            orden->duplicar(),
            elementos,
            inicios_orden_primario,
            indices_orden_primario,
            indices_orden_secundario,
            elemento_nulo,
            filas,
            columnas,
            longitud,
            diccionario.size()
            );
}

template <class T>
void MatrizRala<T>::Iterador::IniciarFilas()
{
    matriz->orden->IniciarFilas(this);
}

template <class T>
void MatrizRala<T>::Iterador::IniciarColumnas()
{
    matriz->orden->IniciarColumnas(this);
}

template <class T>
void MatrizRala<T>::Iterador::SiguienteFila()
{
    matriz->orden->SiguienteFila(this);
}

template <class T>
void MatrizRala<T>::Iterador::SiguienteColumna()
{
    matriz->orden->SiguienteColumna(this);
}

template <class T>
bool MatrizRala<T>::Iterador::QuedanFilas()
{
    return matriz->orden->QuedanFilas(this);
}

template <class T>
bool MatrizRala<T>::Iterador::QuedanColumnas()
{
    return matriz->orden->QuedanColumnas(this);
}

template <class T>
unsigned int MatrizRala<T>::Iterador::FilaActual()
{
    return matriz->orden->FilaActual(this);
}

template <class T>
unsigned int MatrizRala<T>::Iterador::ColumnaActual()
{
    return matriz->orden->ColumnaActual(this);
}

template <class T>
void MatrizRala<T>::Iterador::IniciarOrdenPrimario()
{
    indice_orden_primario = 0;
}

template <class T>
void MatrizRala<T>::Iterador::IniciarOrdenSecundario()
{
    inicio_orden_primario = matriz->inicios_orden_primario[indice_orden_primario];
    fin_orden_primario = matriz->inicios_orden_primario[indice_orden_primario + 1] - 1;
    indice_orden_secundario = inicio_orden_primario;
}

template <class T>
void MatrizRala<T>::Iterador::SiguienteEnOrdenPrimario()
{
    indice_orden_primario++;
}

template <class T>
void MatrizRala<T>::Iterador::SiguienteEnOrdenSecundario()
{
    indice_orden_secundario++;
}

template <class T>
bool MatrizRala<T>::Iterador::QuedanEnOrdenPrimario()
{
    return indice_orden_primario < matriz->longitud_orden_primario;
}

template <class T>
bool MatrizRala<T>::Iterador::QuedanEnOrdenSecundario()
{
    return indice_orden_secundario <= fin_orden_primario;
}

template <class T>
unsigned int MatrizRala<T>::Iterador::OrdenPrimarioActual()
{
    return matriz->indices_orden_primario[indice_orden_primario];
}

template <class T>
unsigned int MatrizRala<T>::Iterador::OrdenSecundarioActual()
{
    return matriz->indices_orden_secundario[indice_orden_secundario];
}

template <class T>
T MatrizRala<T>::Iterador::ElementoActual()
{
    return matriz->elementos[indice_orden_secundario];
}

#endif	/* _MATRIZRALA_H */

