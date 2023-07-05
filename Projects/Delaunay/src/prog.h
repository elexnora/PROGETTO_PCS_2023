#ifndef PROG_H
#define PROG_H
#include <vector>
namespace dalaunay {

class Lati; // la chiamo cosi Point la legge

class Point{
  public:
 double x;
 double y;
 std::vector<Lati *  > edges={}; // vettore di puntatori a lati


Point(double x, double y) : x(x), y(y) {} //costruttore
Point() : x(0), y(0) {} //costruttore vuoto

//  possono essere utili
 Point operator-(const Point & right)
    {
       return Point(x - right.x,y - right.y);
    }
 Point operator+(const Point & right)
    {
       return Point(x + right.x,y + right.y);
    }


 static double dist(Point& p1,Point&p2) // funzione che calcola la distanza tra due punti
 {
    return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
 }

 Lati * isconnect(Point*p); // funzione che mi dice quale lato collega il punto a p

};



class Lati{
public:
    Point * points[2]; // attributo
    Lati(){} // costruttore vuoto
    Lati(Point* p1,Point* p2) // costruttore
    {
        points[0]=p1;
        points[1]=p2;
    }

};



class Triangle{
  public:
 // attributi
 Point * points[3];
 Lati * sides[3];


 Triangle(Point *p1,Point *p2,Point *p3) // costruttore
 {
    points[0]=p1;
    points[1]=p2;
    points[2]=p3;
    ordinapunti();

 }
 Triangle(){}; // costruttore vuoto

void ordinapunti () // restituisce il triangolo ordinato in senso antiorario
 {
    Point AB=*points[1]-*points[0];
    Point AC=*points[2]-*points[0];
    double z =AB.x*AC.y-AB.y*AC.x;
    if (z<0)
    {       
        Point * tmp= points[2];
        points[2]=points[1];
        points[1]=tmp;
    }

 }


 bool ispointin(const Point &p) // passo il puntatore e quindi potrei modificare l'originale, più veloce la chiamata
 {
    // [2] DA FARE --> vedere suggerimento 5

    if(!isInsideCircumcircle(*points[0],*points[1],*points[2],p))
        return false;
    else
    {

        // trovo i due punti più vicini a b
        Point * a= points[0];
        Point * b= points[1];
        //dentro (  *a, * b,  p);
        if(! dentro( *a, * b,  p))
        {
            return false;
        }
        else
            return true;
    }


 }


bool isInsideCircumcircle(const Point& q, const Point& a, const Point& b, const Point& c);


bool dentro ( Point a, Point b, Point c) // restituisce il triangolo ordinato in senso antiorario
 {
    Point AB=a-b;
    Point AC=a-c;
    double z =AB.x*AC.y-AB.y*AC.x;
    if (z<0)
        return false;
    else
        return true;
 }



};
void Dalaunay();

}




#endif // PROG_H

