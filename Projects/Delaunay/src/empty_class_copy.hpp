#ifndef __EMPTY_H
#define __EMPTY_H

#include <algorithm>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>

using namespace std ;


namespace dalaunay{

class Lati; // la chiamo cosi Point la legge

class Point{
public:
 double x;
 double y;
 vector<Lati *  > edges={}; // vettore di puntatori a lati di cui il punto fa parte


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
 inline bool operator==( const Point p2) const
 {
    return (x==p2.x && y==p2.y);
 }
 inline bool operator!=( const Point p2) const
 {
    return (x!=p2.x || y!=p2.y);
 }

 static double dist(Point& p1,Point&p2) // funzione che calcola la distanza tra due punti
 {
    return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
 }

 Lati * isconnect(Point*p); // funzione che mi dice quale lato collega il punto a p

};


class Triangle;

class Lati{
public:
    Point * points[2]; // attributo
    vector<Triangle *  > triangoli={};
    Lati(){} // costruttore vuoto
    Lati(Point* p1,Point* p2) // costruttore
    {
        points[0]=p1;
        points[1]=p2;
        // aggiungo il nuovo lato alla lista di lati di cui il punto fa parte
        points[0]->edges.push_back(this);
        points[1]->edges.push_back(this);

    }
    Lati(Point* p1,Point* p2,int i) // costruttore lato di prova
    {
        points[0]=p1;
        points[1]=p2;

    }


    inline bool operator==( const Lati l) const
    {
       return ((points[0]==l.points[0] || points[0]==l.points[1])||points[1]==l.points[0] || points[1]==l.points[1]);
    }
};



class Triangle{
  public:
 // attributi
 Point * points[3];
 Lati * sides[3];
 vector<Triangle *  > triangoli_adiacenti={};

 Triangle(Point *p1,Point *p2,Point *p3) // costruttore
 {
    points[0]=p1;
    points[1]=p2;
    points[2]=p3;
    ordinapunti();
    Lati* l1;
    Lati* l2 ;
    Lati* l3 ;
    l1= new Lati(points[0],points[1]);
    l2= new Lati(points[1],points[2]);
    l3= new Lati(points[2],points[0]);

    sides[0]=l1;
    sides[1]=l2;
    sides[2]=l3;
    //delete l1;
    //delete l2;
    //delete l3;


 }
 Triangle(Lati *l1,Lati *l2,Lati *l3) // costruttore
 {
    sides[0]=l1;
    sides[1]=l2;
    sides[2]=l3;
    l1->triangoli.push_back(this);
    l2->triangoli.push_back(this);
    l3->triangoli.push_back(this);
    points[0]=l1->points[0];
    points[1]=l2->points[0];
    points[2]=l3->points[0];
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

 inline bool operator==( const Triangle t2) const
 {
    return ((points[0]==t2.points[0]||points[0]==t2.points[1]||points[0]==t2.points[2])&& (points[1]==t2.points[0]||points[1]==t2.points[1]||points[1]==t2.points[2])&& (points[2]==t2.points[0]||points[2]==t2.points[1]||points[2]==t2.points[2]));
 }

 };




 void InputFile(string inputFilePath, vector<Point> &points);
 void Distanzamax(vector<Point> points,Point& p1,Point& p2,Point& p3);
 bool dentro ( Point a, Point b, Point c);
 bool isInsideCircumcircle(Triangle t, Point& q);
 bool ispointin(Triangle t, Point &p);
 bool siintersecano( Point p1, Point p2, Point p3, Point p4);
 void creotriangoli(vector<Triangle> &triangles,int j, Point* p,vector<Lati> &lati);
// void creolati(Triangle &t);
 void Flip(Triangle& ABC, Triangle& CBD, Lati& oldlato);
 //bool isDelaunayValid( Point& a,  Point& b,  Point& c,  Point& d);
 void Posizione(vector<Point*>& point_da_aggiungere, vector<Triangle>& triangles,vector<Lati*>& lati_di_bordo, vector<Lati>& lati, vector<Point*>& punti_di_bordo);
 void ipotesi(Point* punto,vector<Lati>&lati);
 void output (vector<Lati>lati);
//void Delaunay(string InputFilePath);
 void Delaunay();
}



#endif //__EMPTY_H
