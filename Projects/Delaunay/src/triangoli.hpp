#ifndef __TRIANGOLI_H
#define __TRIANGOLI_H


#include <algorithm>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

#include <cmath>
#include <vector>
#include <map>
#include <Eigen/Eigen>
#include<list>
#include <iostream>
#include "empty_class.hpp"


using namespace std ;
using namespace Eigen;



namespace dalaunay
{
bool sinistra (const Point& a, const Point& b, const Point& c);
bool dentro_cerchio_circoscritto(const Triangle& t, const Point& q);
bool dentro_triangolo(const Triangle &t,const Point &p);
void nuovi_triangoli(Triangle* &t, const Point &p,vector <Lati> &lati,vector<Triangle*> &triangoli,int& t_id,int& l_id);
Triangle* trova_triangolo(const Point &p,const vector<Triangle*> &triangoli );
bool verifica_adiacenza(Triangle* &t1,Triangle* &t2,Lati &l);
void rimuovi_lato(vector<Lati>& vec,const Lati &val);
void rimuovi_triangolo(vector<Triangle*>& vec, Triangle* &ptr);
}


#endif //__TRIANGOLI_H
