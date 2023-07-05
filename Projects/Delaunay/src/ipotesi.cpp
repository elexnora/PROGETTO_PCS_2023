#include "ipotesi.hpp"
#include "empty_class.hpp"
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

void nuovo_lato_flip(const Triangle &t1, const Triangle &t2, const Lati &latovecchio,Lati & latonuovo,int& l_id)  // trova il nuovo lato che serve per il flip
{
   Point a;
   Point b;
   for(int i=0;i<3;i++)
   {
       if(t1.points[i]!=latovecchio.points[0] && t1.points[i]!=latovecchio.points[1]) //se il punto del triangolo 1 è diverso da quelli del lato vecchio
           a = t1.points[i];
       if(t2.points[i]!=latovecchio.points[0] && t2.points[i]!=latovecchio.points[1])  //se il punto del triangolo 2 è diverso da quelli del lato vecchio
           b = t2.points[i];
   }
   latonuovo = Lati(a,b,l_id);  // creo il lato con i due punti trovati
   l_id++;  // incremento id perchè è una funzione che chiamo più volte e quindi deve aggiornare l'id del lato nuovo creato
}


void eredita_adiacenze_flip(Triangle* &old,Triangle* &nuovo)  //fa ereditare le adiacenze
{
    Lati lato;  // creo un lato
    if(verifica_adiacenza(old,nuovo,lato))  // se i due triangoli sono adiacenti
      {
       int a;
       for(int i=0;i<3;i++)
        {
          for (int j=0;j<3;j++)
          {
            if(nuovo->sides[j]==old->sides[i] && nuovo->sides[j]==lato)  // se un lato del nuovo triangolo è uguale a uno del vecchio ed è uguale a quello in comune
             {
                nuovo->ad[j] = old->ad[i];  // inserisco quel lato negli adiacenti a quelli del  nuovo triangolo
                a = i;  // mi salvo la posizione
             }
           }
         }

       // viceversa
       for (int j=0; j<3; j++)
           {
           if(old->ad[a]!=nullptr)
           {
             if(old->ad[a]->ad[j]==old)
              {
                old->ad[a]->ad[j] = nuovo;
              }
           }
           }
         }
}


void Flip(Triangle* &ABC,Triangle* &CBD,const Lati &latovecchio,vector<Triangle*>& triangoli,vector<Lati>&lati,int& t_id,int& l_id,const Point &p)
{  // devo trovare il nuovo lato
   rimuovi_lato( lati, latovecchio);  // se faccio il flip devo rimuovere latovecchio dalla lista lati
   Lati latonuovo;
   nuovo_lato_flip(*ABC,*CBD,latovecchio,latonuovo,l_id);  // trova il nuovo lato che si deve formare
   lati.push_back(latonuovo);  //aggiungo alla lista lati il nuovo lato creato

   // creo nuovi triangoli
     Triangle* ADB = new Triangle(latonuovo.points[0], latonuovo.points[1], latovecchio.points[0],t_id,l_id );  // creo il nuovo triangolo con i punti del nuovo lato e uno del vecchio lato
     t_id++;l_id++;l_id++;l_id++;  // incremento il numero dei lati di 3 e il numero del triangolo di uno
     Triangle* ADC = new Triangle(latonuovo.points[0], latonuovo.points[1], latovecchio.points[1],t_id,l_id );
     t_id++;l_id++;l_id++;l_id++;
     triangoli.push_back(ADB);  // inserisco i due nuovi triangoli nella lista
     triangoli.push_back(ADC);

   // aggiorno le adiacenze

     eredita_adiacenze_flip(ABC,ADB);
     eredita_adiacenze_flip(ABC,ADC);
     eredita_adiacenze_flip(CBD,ADB);
     eredita_adiacenze_flip(CBD,ADC);

     rimuovi_triangolo(triangoli, ABC);  // rimuovo i triangoli che c'erano prima del flip
     rimuovi_triangolo(triangoli, CBD);
   // creo le adiacenze tra i due nuovi triangoli creati

     for (int i=0;i<3;i++)
     {
         if(ADB->sides[i]==latonuovo) // se il lato i del nuovo triangolo ABD è quello in comune
             ADB->ad[i] = ADC;  // metto come adiacente il triangolo ADC
         if(ADC->sides[i]==latonuovo)
             ADC->ad[i] = ADB;
     }

     verifica_ipotesi_Delaunay(p,ADB,triangoli,lati,t_id,l_id);  // verifivo che sia rispettata l'ipotesi di Delaunay per ogni triangolo aggiunto
     verifica_ipotesi_Delaunay(p,ADC,triangoli,lati,t_id,l_id);
 }



void verifica_ipotesi_Delaunay(const Point &punto,Triangle* &t,vector<Triangle*> &triangoli,vector<Lati> &lati,int& t_id,int& l_id)
{
   vector<Triangle*> triangoli_da_verificare;
   vector<int> posizione_lato;
   Triangle* T;
   for(int i=0;i<3;i++)
   {
       if(t->ad[i]!=nullptr)  // se ha dei triangoli adiacenti
       {
           if(t->ad[i]->points[0]!=punto && t->ad[i]->points[1]!=punto && t->ad[i]->points[2]!=punto)  // se il punto non fa parte del triangolo in posizione adiacente
            {
               triangoli_da_verificare.push_back(t->ad[i]);  // inserisco nei triangoli da verificare il triangolo i adiacente a t
               posizione_lato.push_back(i);
            }
       }
   }
   while (triangoli_da_verificare.size()>0)  // finchè ci sono triangoli da verificare
   {
        T = triangoli_da_verificare.back();  // parto dal fondo
        Lati l= t->sides[posizione_lato.back()];
        if(dentro_cerchio_circoscritto(*T,punto)) // verifica l'ipotesi di Delaunay
        {
            Flip(t,T,l,triangoli,lati,t_id,l_id,punto);
        }
        triangoli_da_verificare.pop_back(); // dopo averlo verificato lo elimino dalla lista
        posizione_lato.pop_back();
   }
}

}
