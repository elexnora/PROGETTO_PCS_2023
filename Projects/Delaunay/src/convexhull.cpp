#include "empty_class.hpp"
#include "convexhull.hpp"
#include "ipotesi.hpp"
#include "triangoli.hpp"
#include "sorting.hpp"

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


namespace dalaunay {

bool UpperLine(const Point& p1,
               const Point& p2,
               const Point& p3)
{
  double dx1 = p2.x - p1.x;
  double dy1 = p2.y - p1.y;
  double dx2 = p3.x - p1.x;
  double dy2 = p3.y - p1.y;
  if (dx1*dy2 > dy1*dx2)
    return true;
  else if (dx1*dy2 < dy1*dx2)
    return false;
  else
    return (dx1*dx2+dy1*dy2 >= TOLL);
}




   vector<Point> ConvexHull(vector<Point> &points)
  {
    if(points.size() < 2)
       return points;
    else
    {

    SortLibrary::MergeSort(points,0,int(points.size()-1));

    // eseguo il lowerHull
    vector<Point> lowerHull;
    for(unsigned int i=0;i<points.size();++i)
    {
      while(lowerHull.size()>= 2 && !UpperLine(lowerHull[lowerHull.size()-2],lowerHull[lowerHull.size()-1],points[i]))
        {
          lowerHull.pop_back();
        }
      lowerHull.push_back(points[i]);
    }

    // eseguo l'upperHull
    vector<Point> upperHull;
    for(int i=points.size()-1;i>=0;--i)
    {
      while(upperHull.size() >= 2 && !UpperLine(upperHull[upperHull.size()-2],upperHull[upperHull.size()-1],points[i]))
        {
          upperHull.pop_back();
        }
          upperHull.push_back(points[i]);
    }

    // rimuovo il primo e l'ultimo punto dall'upperHull
    upperHull.pop_back();
    upperHull.erase(upperHull.begin());

    // concateno i due hulls
    vector<Point> hull(lowerHull);
    hull.insert(hull.end(), upperHull.begin(), upperHull.end());

    return hull;
  }
 }


    vector<Triangle*> TriangoliConvexHull (const vector<Point>& punti,vector <Lati>& lati,int& t_id,int& l_id)
    {
         vector<Triangle*> triangoli;
         int size = int(punti.size());
         for(int i=0;i<size;i++)
         {
           Lati l = Lati(punti[i],punti[(i+1)%size],l_id);
           l_id++;
           lati.push_back(l);
         }
         for(int i=1;i<size-1;i++)
         {
             Triangle* nuovoTriangolo =new Triangle(punti[0],punti[i],punti[(i+1)],t_id,l_id);
             t_id++;l_id++;l_id++;l_id++;

             if(i==1)
             {
                triangoli.push_back(nuovoTriangolo);
             }

             else
             {
             Lati l = Lati(punti[0],punti[i],l_id);
             l_id++;
             lati.push_back(l);
             triangoli.push_back(nuovoTriangolo);
             for(int j=0;j<int(triangoli.size()-1);j++)
             {                
               if(verifica_adiacenza(nuovoTriangolo,triangoli[j],l))
               {
                   for(int k=0;k<3;k++)
                   {
                     if(nuovoTriangolo->sides[k]==l)
                      {
                       nuovoTriangolo->ad[k]=triangoli[j];
                       }
                     if(triangoli[j]->sides[k]==l)
                       {
                        triangoli[j]->ad[k]=nuovoTriangolo;                       
                       }
                   }
               }
            }
            verifica_ipotesi_Delaunay(nuovoTriangolo->points[2],nuovoTriangolo,triangoli,lati,t_id,l_id);
            }
         }
       return triangoli;
}

}
