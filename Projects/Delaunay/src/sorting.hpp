#ifndef __SORTING_H
#define __SORTING_H

#include "iostream"
#include "list"
#include "Eigen/Eigen"
#include "map"
#include <vector>

using namespace std;
using namespace Eigen;

namespace SortLibrary {

  template<typename T>  // le function template vengono definite nell'hpp
  void Merge(vector<T>& v,
             const unsigned int& sx,
             const unsigned int& cx,
             const unsigned int& dx)
  { unsigned int i = sx; unsigned int j = cx+1; unsigned int k = 0;

    vector<T> b(v.size());

    while ((i <= cx) && (j <= dx))
    {
      if (v[i] <= v[j])
      {
        b[k] = v[i]; i = i+1;
      }
      else
      {
        b[k] = v[j]; j = j+1;
      }
      k = k+1;
      }
      for ( ; i <= cx; i = i+1, k = k+1)
          b[k] = v[i];
      for ( ; j <= dx; j = j+1, k = k+1)
          b[k] = v[j];
      for (i = sx; i <= dx; i = i+1)
          v[i] = b[i-sx];



  }

  template<typename T>
  void MergeSort(vector<T>& v,
                 const unsigned int& sx,
                 const unsigned int& dx)

  {if (sx < dx)
     {
     unsigned int centro = (sx+dx)/2;
     MergeSort( v, sx, centro );
     MergeSort( v, centro+1, dx );
     Merge( v, sx, centro, dx );
     }

  }

}

#endif // __SORTING_H
