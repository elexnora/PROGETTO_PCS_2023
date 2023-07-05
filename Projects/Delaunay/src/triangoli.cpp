#include "empty_class.hpp"
#include "triangoli.hpp"
#include "ipotesi.hpp"

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

bool sinistra (const Point& a,const Point& b,const Point &c)   // ritorna falso se z è minore di zero
{   // mi dice se c sta a sinistra del segmento ab

    Point B=b;// creo una copia dei punti non const per poter usare l'operator
    Point C=c;
    Point AB = B-a;
    Point AC = C-a;
    double z = AB.x*AC.y-AB.y*AC.x;
    if (z< 0)
        return false;
    else
        return true;
}


 bool dentro_triangolo(const Triangle &t,const Point &p)   // dice se il punto p sta dentro il triangolo t
 {
    if(!dentro_cerchio_circoscritto(t,p))
        return false;
    else
    {
        Point a = t.points[0];
        Point b = t.points[1];
        Point c = t.points[2];
        if(sinistra(a,b,p) && sinistra(b,c,p) && sinistra(c,a,p))   // se il punto sta a sinistra di tutti i lati
        {
            return true;
        }
        else
            return false;
    }
 }


 bool dentro_cerchio_circoscritto(const Triangle &t, const Point &q)
 {//mi dice se il punto p sta dentro il cerchio circoscritto del triangolo t
     Point a = t.points[0];
     Point b = t.points[1];
     Point c = t.points[2];

     double ax = a.x - q.x;
     double ay = a.y - q.y;
     double bx = b.x - q.x;
     double by = b.y - q.y;
     double cx = c.x - q.x;
     double cy = c.y - q.y;

     double determinante = ax * (by * (cx*cx+cy*cy) - cy*(bx*bx+by*by)) - bx * (ay*(cx*cx+cy*cy) - cy*(ax*ax+ay*ay)) + cx * (ay*(bx*bx+by*by) - by*(ax*ax+ay*ay));

     return determinante > 0;

 }


void nuovi_triangoli(Triangle* &t,const Point& p,vector <Lati> &lati,vector<Triangle*> &triangoli,int& t_id,int& l_id)  // crea i 3 triangoli interni al triangolo t
{
    Triangle old=*t;   //devo eliminare t da triangoli
    Triangle* t1=new Triangle(old.sides[0].points[0],old.sides[0].points[1],p,t_id,l_id);
    triangoli.push_back(t1);  //aggiunge il triangolo formato dal punto p e dai punti del primo lato del triangolo vecchio
    t_id++;l_id++;l_id++;l_id++;  // incremento id dei triangoli di 1 e dei lati di 3
    Triangle* t2=new Triangle(old.sides[1].points[0],old.sides[1].points[1],p,t_id,l_id);
    triangoli.push_back(t2);
    t_id++;l_id++;l_id++;l_id++;
    Triangle* t3=new Triangle(old.sides[2].points[0],old.sides[2].points[1],p,t_id,l_id);
    triangoli.push_back(t3);
    t_id++;l_id++;l_id++;l_id++;
    Triangle* tr[3]={t1,t2,t3};
    for(int i=0;i<3;i++)
       {
           if(t->ad[i]!=nullptr)  // se il vecchio triangolo ha un triangolo adiacente
           {
               Lati l;
               if(verifica_adiacenza(t,tr[i],l))  // verifico adiancenze tra t old e i nuovi triangoli, l lato di adiacenza
               {

               for(int j=0;j<3;j++)
               {
                   if( tr[i]->sides[j]==l)
                   {
                        tr[i]->ad[j] = t->ad[i];
                   }

                   if(old.ad[i]->ad[j]==t)
                   {
                      old.ad[i]->ad[j] = tr[i];
                   }    // aggiorno l'adiacenza dell'altro

               }
               }
           }

           lati.push_back(Lati(old.points[i],p,l_id));
           l_id++;
       }
       rimuovi_triangolo(triangoli,t);

       for(int i=0;i<3;i++)
       {
          for (int j=0;j<3;j++)
          {
              if(t1->sides[j]==t2->sides[i])
              {
                  t1->ad[j] =t2;
                  t2->ad[i] = t1;
              }
              if(t1->sides[j]==t3->sides[i])
              {
                  t1->ad[j] =t3;
                  t3->ad[i] = t1;
              }
              if(t2->sides[j]==t3->sides[i])
              {
                  t2->ad[j] =t3;
                  t3->ad[i] = t2;
              }

          }
       }

       verifica_ipotesi_Delaunay(p,t1,triangoli,lati,t_id,l_id);
       verifica_ipotesi_Delaunay(p,t2,triangoli,lati,t_id,l_id);
       verifica_ipotesi_Delaunay(p,t3,triangoli,lati,t_id,l_id);
   }

  Triangle* trova_triangolo(const Point &p,const vector<Triangle*> &triangoli )  // dice l'id del triangolo in cui il punto p è dentro
  {

    for (int j=0;j<int(triangoli.size());j++)
    {
        if(dentro_triangolo(*triangoli[j],p))
        {
            return triangoli[j];
        }
    }
    return nullptr;
  }

 bool verifica_adiacenza(Triangle* &t1, Triangle* &t2,Lati &l)  //dati due triangoli mi dice se sono adiacenti e mi restituisce il lato tra essi
 {
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            if(t1->sides[i]==t2->sides[j])  // cicla sui lati dei due triangoli, se sono uguali restituisce true e mi salva il lato
            {
                l = t1->sides[i];
                return true;
            }
        }
    }
    return false;
 }

 void rimuovi_lato(vector<Lati>& vec,const Lati &val)// rimuove l'elemento val dalla lista vec
 {
     for (auto it=vec.begin();it!=vec.end();++it)
     {
         if (*it==val)
         {
             vec.erase(it);
             break;
         }

     }
 }

 void rimuovi_triangolo(vector<Triangle*>& vec, Triangle* &ptr)  // rimuove l'elemento val dalla lista vec
 {
     auto it = std::find(vec.begin(), vec.end(), ptr);
     if (it != vec.end())
     {
         vec.erase(it);
         delete ptr;
     }

 }

}
