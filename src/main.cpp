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

#define INDICE_ARG_METODO  1
#define INDICE_ARG_PAGINAS 2

using namespace std;

int main(int argc, char ** argv)
{
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
        cerr << "Error: No se pudo abrir el archivo de links: '"<< argv[INDICE_ARG_PAGINAS] << "'." << endl;
        return 1;
    }

    getline(archivo_entrada, linea);
    c_linea = new char[linea.size()];
    strcpy(c_linea, linea.c_str());
    unsigned dimensiones = atoi(c_linea);
    // FIXME No se como ir a la siguiente linea
    // para obtener la cantidad de links
    unsigned links = atoi(c_linea);
    delete c_linea;
    archivo_entrada.close();

    cout << "Cantidad de sitios: " << dimensiones << endl;
    cout << "Cantidad de links: " << links << endl;

    //ordenes de matrices generales
    MatrizRala<float>::Orden* orden_por_filas = new MatrizRala<float>::OrdenPorFilas();
    MatrizRala<float>::Orden* orden_por_columnas = new MatrizRala<float>::OrdenPorColumnas();

    float cero = 0.0f;
    float uno = 1.0f;
    // creador de la matriz W, valor por default 0
    MatrizRala<float>::Creador creador_W(orden_por_filas, dimensiones, dimensiones, cero);

    // Se llena la matriz
    char* valor;
    int i, j;
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
        creador_W.Agregar(j - 1, i - 1, uno);
    }

    archivo_entrada.close();

    // Construi la Matriz W
    const MatrizRala<float>* W = creador_W.Crear();

    const MatrizRala<float>* W_por_filas = W->OrdenadaPorFilas();
    typename MatrizRala<float>*::Iterador iterador(W_por_filas);

    // Contruimos P'
    // En un paso intermedio obtendremos P
    // primero obtengo los grados iterando W
    vector<float> grados(dimensiones);
    int index = 0;
    for (iterador.IniciarFilas(); iterador.QuedanFilas(); iterador.SiguienteFila())
    {
        grados[iterador.FilaActual()] = 0;
        for (iterador.IniciarColumnas(); iterador.QuedanColumnas(); iterador.SiguienteColumna())
        {
            grados[iterador.FilaActual()] += iterador.ElementoActual();
        }
    }
    delete W_por_filas;

    // Ahora armamos P' iterando nuevamente W
    typename MatrizRala<float>*::Iterador iterador2(W_por_filas);

    MatrizRala<float>::Creador creador_P1(orden_por_filas, dimensiones, dimensiones, cero);
    for (iterador2.IniciarFilas(); iterador2.QuedanFilas(); iterador2.SiguienteFila())
    {
        for (iterador2.IniciarColumnas(); iterador2.QuedanColumnas(); iterador2.SiguienteColumna())
        {
            creador_P1.Agregar(iterador2.ColunmaActual(), iterador2.FilaActual(), uno / grados[iterador2.ColunmaActual()])
        }
    }
    delete W_por_filas;

    // Tenemos ahora que P1 == P
    // para que P1 = P + D, pondremos en
    // la diagonal 1/n
    for (int i = 0; i < dimensiones; ++i)
        creador_P1.Agregar(i, i, uno / dimensiones);

    const MatrizRala<float>* P1 = creador_P1.Crear();

    return 0;
}
