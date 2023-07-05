#include "empty_class.hpp"
#include "convexhull.hpp"
#include "funzioni.hpp"
#include "triangoli.hpp"

#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>
#include <unordered_set>




using namespace std;

namespace dalaunay
{
void Delaunay()
{

    string PercorsoFileInput = "C:/Users/matteo/OneDrive/Desktop/PCS2023_Exercises/Projects/Delaunay/Dataset/Test2.csv";
    vector<Point> punti;

    // leggo i dati dal file
    LeggoFileInput(PercorsoFileInput,punti);

    // prendo il convexhull
    int id_triangoli = 0;
    int id_lati = 0;

    vector<Triangle*> triangoli;
    vector <Lati> lati; // creo un vettore di lati
    vector<Point> copia_punti=punti;
    vector <Point> convexhull =  ConvexHull(copia_punti);  // dati i punti mi restituisce il Convehull di tali punti
    triangoli = TriangoliConvexHull(convexhull,lati,id_triangoli,id_lati);

    vector<Point> punti_da_aggiungere;
    for (int i=0;i<int(punti.size());i++)
    {
        bool a= true;
         for (int j=0;j<int(convexhull.size());j++)
         {
            if(punti[i]==convexhull[j])
            {
                a=false;
                break;
            }
         }
         if(a)
            {punti_da_aggiungere.push_back(punti[i]);}
    }


    // cerco dentro a quale triangolo si trova il punto
    for(int i=0;i<int(punti_da_aggiungere.size());i++)
    {
        Triangle* T=trova_triangolo(punti_da_aggiungere[i],triangoli);
        nuovi_triangoli(T,punti_da_aggiungere[i],lati,triangoli,id_triangoli,id_lati);
    }
    // output schermo
    cout<<"Lati: "<<endl;
    for (int i=0;i<int(lati.size());i++)
    {
        Point p1;
        p1 = lati[i].points[0];
        Point p2;
        p2 = lati[i].points[1];
        cout<<to_string(i+1)<<") "<<"punto di inizio: ("<<p1.x<<","<<p1.y<<")       punto di fine: ("<<p2.x<<","<<p2.y <<")"<<endl;
    }
    // output file
    string output = "C:/Users/matteo/OneDrive/Desktop/PCS2023_Exercises/Projects/Delaunay/Dataset/output.csv";
    CreoFileOutput(output,lati);

    for (size_t i = 0; i < triangoli.size(); ++i)
    {
        delete triangoli[i];  // Dealloca la memoria puntata da ogni puntatore
    }
    triangoli.clear();
}

}
