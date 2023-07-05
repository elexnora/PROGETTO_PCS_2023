#include "prog.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

#include <cmath>


using namespace std;
using namespace dalaunay;

vector<Lati> sides;  //creo una lista di lati

void Dalaunay()
{
    // LEGGO I DATI


    ifstream file;
    file.open("./test1.csv");
    string line;
    vector<Point> points; //creo una lista di punti
    getline(file,line);
    while (getline(file, line))//while file >> file
    {
       stringstream ss; //converter
       ss << line;
       int id;
       ss>>id;
       double x,y;
       ss>>x>>y;
       points.push_back(Point(x,y)); // creo le istanze Punto nella lista points
    }
    file.close();


    // PRENDO IL TRIANGOLO DI AREA MAGGIORE


    //mergesort ??? o meglio trovare il triangolo di area meaasima O(n^3)????

    // PRENDO I DUE PUNTI PIù LONTANI

    //confronto tutte le coppie O(n^2)
    double max;
    Point * p1;
    Point * p2;
    max=0;
    for (int i=0;i<points.size();i++)
    {
        for (int j=0;j<points.size();j++)
        {
            if (j!=i)
            {
               double distance=Point::dist(points[i],points[j]);
               if(distance>max)
               {
                   max=distance;
                   p1=&points[i];
                   p2=&points[j];

               }
            }
        }

    }


    // [1] DA FARE:
    //devo trovare il punto p3 alla massima disanza punto retta per p1 p2


    Point * p3;


    vector<Triangle> triangles;  // creo una lista di triangoli
    triangles.push_back(Triangle(p1,p2,p3)); // aggiungo il primo triangolo

    vector<Point*> punti_di_bordo; // andrebbe poi aggiornata!!!!!
    punti_di_bordo.push_back(p1);
    punti_di_bordo.push_back(p2);
    punti_di_bordo.push_back(p3);

    //creo una lista di punti senza p1 p2 p3
    vector<Point*> point_da_aggiungere;
    for (int i=0;i<points.size();i++)
    {
        if(&points[i]==p1||&points[i]==p2||&points[i]==p3)
        {
            continue;
        }
        else
        {
            point_da_aggiungere.push_back(&points[i]);
        }
    }


    // CONTROLLO LA POSIZIONE DEL NUOVO PUNTO RISPETTO AI TRIANGOLI


    // [2] DA FARE: ( in hpp ispointin ) --> vedere suggerimento 5





    for (int i=0;i<point_da_aggiungere.size();i++)// per ogni punto
    {
        bool in=false;

        for (int j=0;j<triangles.size();j++) // per ogni triangolo
        {
            if(triangles[j].ispointin(*point_da_aggiungere[i])) // se il punto sta dentro il triangolo (dereferenzio (*))
            {
                // creo 3 nuovi triangoli
                in=true;
                Triangle old=triangles[j];
                triangles[j]=Triangle(old.points[0],old.points[1],point_da_aggiungere[i]);
                triangles.push_back(Triangle(old.points[1],old.points[2],point_da_aggiungere[i]));
                triangles.push_back(Triangle(old.points[2],old.points[0],point_da_aggiungere[i]));

                // devo verificare ipotesi dalaunay?????

                break; // esco dal ciclo for
            }

        }


        if(!in) // se il punto è esterno a tutti
        {
           for (int j=0;j<triangles.size();j++)  // per ogni triangolo
           {
               for(int k=0;k<3;j++)
               {
                Lati l= Lati(point_da_aggiungere[i],triangles[j].points[k]);


                // [3] DA FARE : ( sotto in creolati )
                // sfutto la lista punti_di_bordo
                //verifico che non crea intersezioni
                // se non interseca aggiungo


                sides.push_back(l);

               }

               //fare una funzione crea triangoli al posto di tre if????

               if ((point_da_aggiungere[i]->isconnect(triangles[j].points[0]))&&(point_da_aggiungere[i]->isconnect(triangles[j].points[1])))
               {
                   triangles.push_back(Triangle(point_da_aggiungere[i],triangles[j].points[0],triangles[j].points[1]));
               }
               else if ((point_da_aggiungere[i]->isconnect(triangles[j].points[0]))&&(point_da_aggiungere[i]->isconnect(triangles[j].points[2])))
               {
                   triangles.push_back(Triangle(point_da_aggiungere[i],triangles[j].points[0],triangles[j].points[2]));
               }
               else if ((point_da_aggiungere[i]->isconnect(triangles[j].points[1]))&&(point_da_aggiungere[i]->isconnect(triangles[j].points[2])))
               {
                   triangles.push_back(Triangle(point_da_aggiungere[i],triangles[j].points[1],triangles[j].points[2]));
               }
           }
        punti_di_bordo.push_back(point_da_aggiungere[i]);
        }
    }





}




dalaunay::Lati * Point::isconnect(Point*p) // verifico se è collegato con p
{
    for (auto l : edges)// per ogni lato l
    {
       if((l->points[0]==this&&p==l->points[1])||(l->points[1]==this&&p==l->points[0]))
       {
           return l;
       }

    }
    return nullptr;
}








void creolati(Triangle &t)
{
    Lati * l1=t.points[0]->isconnect(t.points[1]);
    if (l1==nullptr) // se non eseiste il lato
    {

        sides.push_back(Lati(t.points[0],t.points[1])); // creo i nuovi lati del triangolo


        // [3] DA FARE:
        //creare funzione che verifica le intersezioni
        // aggiungo il lato e il nuovo triangolo


    }


}



//////////////////////////////////////////////////////////////////////////////


// DA VERIFICARE

// Verifica se un punto Q è all'interno del circumcerchio di un triangolo ABC
bool isInsideCircumcircle(const Point& q, const Point& a, const Point& b, const Point& c)
{

    double ax = a.x - q.x;

    double ay = a.y - q.y;

    double bx = b.x - q.x;

    double by = b.y - q.y;

    double cx = c.x - q.x;

    double cy = c.y - q.y;



    double det = ax * (by * (cx*cx+cy*cy) - cy*(bx*bx+by*by)) - bx * (ay*(cx*cx+cy*cy) - cy*(ax*ax+ay*ay)) + cx * (ay*(bx*bx+by*by) - by*(ax*ax+ay*ay));
    // dovrebbe essere giusta


    return det > 0.0;

}


// ???

// Verifica l'ipotesi di Delaunay tra due triangoli adiacenti

bool isDelaunayValid(const Point& a, const Point& b, const Point& c, const Point& d) {

    // Calcola i quadrati delle distanze tra i punti

    double abSquared = std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2);

    double acSquared = std::pow(a.x - c.x, 2) + std::pow(a.y - c.y, 2);

    double adSquared = std::pow(a.x - d.x, 2) + std::pow(a.y - d.y, 2);

    double bcSquared = std::pow(b.x - c.x, 2) + std::pow(b.y - c.y, 2);

    double bdSquared = std::pow(b.x - d.x, 2) + std::pow(b.y - d.y, 2);



    // Verifica l'ipotesi di Delaunay

    if ((abSquared * acSquared * adSquared) + (abSquared * bcSquared * bdSquared) >

        (acSquared * bcSquared * bdSquared) + (abSquared * acSquared * bdSquared)) {

        return false;

    }



    return true;

}
