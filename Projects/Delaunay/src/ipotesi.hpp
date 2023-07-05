#ifndef __IPOTESI_H
#define __IPOTESI_H

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
void nuovo_lato_flip(const Triangle &ABC,const Triangle &CBD,const Lati &latovecchio,Lati & latonuovo,int& l_id);
void eredita_adiacenze_flip(Triangle*& old,Triangle*& nuovo);
void Flip(Triangle* &ABC,Triangle* &CBD,const Lati& latovecchio,vector<Triangle*>& triangoli,vector<Lati>&lati,int& t_id,int& l_id,const Point &p);
void verifica_ipotesi_Delaunay(const Point &punto,Triangle*& t,vector<Triangle*> &triangoli,vector<Lati> &lati,int& t_id,int& l_id);
}


#endif // __IPOTESI_H
