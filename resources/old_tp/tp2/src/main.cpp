#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stdexcept>
#include <string>
#include <string.h>

#include "MatrizRala.h"
#include "MetodoDirecto.h"
#include "MetodoIterativo.h"

#define INDICE_ARG_METODO 1
#define INDICE_ARG_LINKS 2
#define INDICE_ARG_PAGINAS 3
#define INDICE_ARG_RANKING 4

#define INDICE_ARG_ITERATIVO_ITERACIONES 5
#define INDICE_ARG_ITERATIVO_RESIDUO 6

using namespace std;

template <class T>
const MatrizRala<T>* obtenerD(const MatrizRala<T>* A, typename MatrizRala<T>::Orden* orden, T uno)
{
    typename MatrizRala<T>::Creador creador_D(orden, A->Filas(), A->Columnas(), A->ElementoNulo());

    const MatrizRala<T>* A_por_columnas = A->OrdenadaPorColumnas();

    typename MatrizRala<T>::Iterador iterador(A_por_columnas);

    T sumatoria;

    //busco la suma de todos los elementos
    for(iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
    {
        sumatoria = A->ElementoNulo();

        for(iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
        {
            sumatoria+= iterador.ElementoActual();
        }

        if(sumatoria != A->ElementoNulo())
            creador_D.Agregar(iterador.ColumnaActual(), iterador.ColumnaActual(), uno / sumatoria);
    }

    if(A_por_columnas != A)
        delete A_por_columnas;

    return creador_D.Crear();
}

void mostrarAyudaYSalir()
{
    cerr << "Modo de uso: " << endl
            << "a) -directo "
            << "'archivo de links' "
            << "'archivo de paginas' "
            << "'archivo de ranking'"
            << endl
            << "b) -iterativo "
            << "'archivo de links' "
            << "'archivo de paginas' "
            << "'archivo de ranking'"
            << "'iteraciones' "
            << "'residuo'"
            << endl;

    exit(1);
}

int main(int argc, char ** argv)
{
    if(argc < 2)
        mostrarAyudaYSalir();

    bool directo = strcmp(argv[INDICE_ARG_METODO], "-directo") == 0;
    bool iterativo = strcmp(argv[INDICE_ARG_METODO], "-iterativo") == 0;

    if(!directo && !iterativo)
        mostrarAyudaYSalir();

    if((directo && argc != 5) || (iterativo && argc != 7))
        mostrarAyudaYSalir();

    //buffer de uso general
    char* c_linea;
    //linea de uso general
    string linea;

    //archivo de entrada de uso general
    ifstream archivo_entrada;
    //archivo de salida de uso general
    ofstream archivo_salida;

    // Abro el archivo de paginas, para saber cuantas
    // paginas tengo para poder armar la matriz de las dimensiones necesarios
    archivo_entrada.open(argv[INDICE_ARG_PAGINAS]);

    if(!archivo_entrada.is_open())
    {
        cerr << "Error: No se pudo abrir el archivo de paginas: '"<< argv[INDICE_ARG_PAGINAS] << "'." << endl;
        return 1;
    }

    getline(archivo_entrada, linea);
    c_linea = new char[linea.size()];
    strcpy(c_linea, linea.c_str());
    unsigned dimensiones = atoi(c_linea);
    delete c_linea;
    archivo_entrada.close();

    //cout << "Cantidad de sitios: " << dimensiones << endl;

    //Abro el archivo de links para llenar la matriz
    archivo_entrada.open(argv[INDICE_ARG_LINKS]);

    if(!archivo_entrada.is_open())
    {
        cerr << "Error: No se pudo abrir el archivo de links: '"<< argv[INDICE_ARG_LINKS] << "'." << endl;
        return 1;
    }

    //cantidad de links
    getline(archivo_entrada, linea);
    c_linea = new char[linea.size()];
    strcpy(c_linea, linea.c_str());
    unsigned links = atoi(c_linea);
    delete c_linea;

    //ordenes de matrices generales
    MatrizRala<float>::Orden* orden_por_filas = new MatrizRala<float>::OrdenPorFilas();
    MatrizRala<float>::Orden* orden_por_columnas = new MatrizRala<float>::OrdenPorColumnas();

    //elemento nulo
    float cero = 0.0f;
    float uno = 1.0f;

    //creador de la matriz W
    //la ordeno por filas porque va a multiplicar por izquierda
    MatrizRala<float>::Creador creador_w(orden_por_filas, dimensiones, dimensiones, cero);

    //indices temporales
    int i, j;

    char* valor;

    //por cada link
    for(unsigned link = 0; link < links; link++)
    {
        getline(archivo_entrada, linea);
        c_linea = new char[linea.size()];
        strcpy(c_linea, linea.c_str());

        valor = strtok(c_linea, " ");
        i = atoi(valor);
        valor = strtok(NULL, " ");
        j = atoi(valor);

        delete c_linea;

        //Ignoro los autolinks
        if (i == j)
            continue;

        //La pagina j tiene un link a la pagina i
        creador_w.Agregar(j - 1, i - 1, uno);
    }

    archivo_entrada.close();

    // Construi la Matriz W
    const MatrizRala<float>* W = creador_w.Crear();

    // Vamos a construir a partir de W, la matriz D y generar la matriz A tal
    // que sea (I - p W D)
    const MatrizRala<float>* D = obtenerD(W, orden_por_columnas, uno);

    cout.precision(3);

    //cout << "W" << endl << *W << endl;

    //cout << "D" << endl << *D << endl;

    float p = 0.85f;
    //cout << "p = " << p << endl << endl;

    //la ordeno por filas porque va a multiplicar por izquierda
    const MatrizRala<float>* WD = W->MultiplicarPor(D, orden_por_filas);

    //cout << "WD" << endl << *WD << endl;

    //la ordeno por filas porque va a restar con otra ordenada por filas
    const MatrizRala<float>* pWD = WD->MultiplicarPor(p, orden_por_filas);

    //cout << "pWD" << endl << *pWD << endl;

    //la ordeno por filas porque va a restar con otra ordenada por filas
    const MatrizRala<float>* I = identidad(orden_por_filas, dimensiones, dimensiones, cero);

    //cout << "I" << endl << *I << endl;

    //no importa el orden
    const MatrizRala<float>* I_menos_pWD = I->Restarle(pWD, orden_por_filas);

    //cout << "I - pWD" << endl << *I_menos_pWD << endl;

    //creo un vector b lleno de unos
    //lo ordeno por columnas pq es mejor para el metodo Calcular
    MatrizRala<float>::Creador creador_b(orden_por_columnas, dimensiones, 1, cero);

    for (unsigned i = 0; i < dimensiones; ++ i)
        creador_b.Agregar(i, 0, uno);

    const MatrizRala<float>* b = creador_b.Crear();

    //cout << "b" << endl << *b << endl;

    //metodo que va a calcular el ranking
    CalculadorDeRanking<float>* metodo;

    //valor inicial (solo para el iterativo)
    const MatrizRala<float>* x0;

    //si selecciono el metodo directo
    if(directo)
    {
        metodo = new MetodoDirecto<float>();
    }
    //si selecciono el metodo iterativo
    else if(iterativo)
    {
        unsigned iteraciones = atoi(argv[INDICE_ARG_ITERATIVO_ITERACIONES]);
        float residuo = atof(argv[INDICE_ARG_ITERATIVO_RESIDUO]);

        //creo un vector x0 con un solo uno en la primera posicion
        MatrizRala<float>::Creador creador_x0(orden_por_columnas, dimensiones, 1, cero);
        creador_x0.Agregar(0, 0, uno);
        x0 = creador_x0.Crear();

        //cout << "x0" << endl << *x0 << endl;

        metodo = new MetodoIterativo<float>(iteraciones, residuo, uno, x0);
    }

    //resuelvo el problema con el metodo
    const MatrizRala<float>* x = metodo->Calcular(I_menos_pWD, b);

    //cout << "x" << endl << *x << endl;

    // Abro el archivo de ranking para guardarlo
    archivo_salida.open(argv[INDICE_ARG_RANKING]);

    if(!archivo_salida.is_open())
    {
        cerr << "Error: No se pudo abrir el archivo de ranking: '"<< argv[INDICE_ARG_RANKING] << "'." << endl;
        return 1;
    }

    archivo_salida << *x;

    archivo_salida.close();

    delete W;
    delete D;
    delete WD;
    delete pWD;
    delete I;
    delete I_menos_pWD;
    delete b;
    delete x;
    delete orden_por_columnas;
    delete orden_por_filas;
    delete metodo;

    if(iterativo)
    {
        delete x0;
    }

    return 0;
}
