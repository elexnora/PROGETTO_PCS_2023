#ifndef __FUNZIONI_H
#define __FUNZIONI_H

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

void CreoFileOutput(const string &outputFilePath,const vector<Lati>&mylist);
void LeggoFileInput(const string &inputFilePath,vector<Point> &punti);

}

#endif // __FUNZIONI_H
