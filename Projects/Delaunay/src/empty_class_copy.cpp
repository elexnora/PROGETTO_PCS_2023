#include "empty_class_copy.hpp"
//#include "empty_class.hpp"

#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>


using namespace std;


namespace dalaunay
{

vector<Lati> sides;  //creo una lista di lati

void Delaunay()
{

    //string InputFilePath = "C:/Users/matteo/OneDrive/Desktop/PCS2023_Exercises/Projects/Delaunay/Dataset/Test2.csv";
    string InputFilePath = "C:/Users/matteo/OneDrive/Desktop/PCS2023_Exercises/Projects/Delaunay/Dataset/prova.csv";
    vector<Point> points;

    // LEGGO I DATI
    InputFile(InputFilePath, points);

    // PRENDO IL TRIANGOLO DI AREA MAGGIORE
    Point * p1;
    Point * p2;
    Point * p3;
    //////////////////////////////////
    /*
    Point  punto1=points[0];
    Point  punto2=points[1];
    Point  punto3=points[2];
    points.push_back(punto1);
    points.push_back(punto2);
    points.push_back(punto3);
    */
    ///////////////////////////////////////

    int l=points.size();

    points.push_back(points[0]);
    points.push_back(points[1]);
    points.push_back(points[2]);

    Point &a=points[l];
    Point &b=points[l+1];
    Point &c=points[l+2];

    /*
    // DA MODIFICARE
    Point &a=points[0];
    Point &b=points[1];
    Point &c=points[2];
    */

    /////////////////////////////////////////////////////////////////////////////////////// non prende tutti i punti
    /// // sovrascrive i 3 migliori in points[0],1,2
    Distanzamax(points, a,b,c); // [2] CONVEX HULL
    //points.erase(remove(points.begin(),points.end(),points[0]),points.end());
    p1=&a;
    p2=&b;
    p3=&c;
    points.pop_back();
    points.pop_back();
    points.pop_back();

    vector<Triangle> triangles;  // creo una lista di triangoli
    Triangle t=Triangle(p1,p2,p3);
    triangles.push_back(t); // aggiungo il primo triangolo

    vector<Lati> lati;  // creo una lista di lati
    lati.push_back(*t.sides[0]); // aggiungo i primi 3 lati
    lati.push_back(*t.sides[1]);
    lati.push_back(*t.sides[2]);

    vector<Point*> punti_di_bordo; // andrebbe poi aggiornata!!!!!
    punti_di_bordo.push_back(p1);
    punti_di_bordo.push_back(p2);
    punti_di_bordo.push_back(p3);

    vector<Lati*> lati_di_bordo; // andrebbe poi aggiornata!!!!!
    lati_di_bordo.push_back(triangles[0].sides[0]);
    lati_di_bordo.push_back(triangles[0].sides[1]);
    lati_di_bordo.push_back(triangles[0].sides[2]);

    //creo una lista di punti senza p1 p2 p3
    vector<Point*> point_da_aggiungere;


    for (int i=0;i<int(points.size());i++)
    {

        if(points[i]!=*triangles[0].points[0] && points[i]!=*triangles[0].points[1] && points[i]!=*triangles[0].points[2])
        {
            point_da_aggiungere.push_back(&points[i]);
            //continue;
        }

    }

    // CONTROLLO LA POSIZIONE DEL NUOVO PUNTO RISPETTO AI TRIANGOLI

    // FINO A QUI ARRIVA
//    Posizione(point_da_aggiungere,triangles,lati_di_bordo, lati,punti_di_bordo); // PROBLEMA QUI //////////////////////////////////////
    // QUI NO

    //finita la mesh
    // IPOTESI DI DELAUNAY

    //ipotesi(point_da_aggiungere,lati);


    //[6] OUTPUT

    //output (lati);
    /////////////////////////// NON FUNZIONA ERASE
    //elimino i duplicati
    /*
    for(int i=0;i<int(lati.size());i++)
    {
        if (auto it=find(lati.begin(),lati.end(),lati[i])!=lati.end())
        {
            lati.erase(lati.begin()+it);
        }


    }
    */



    ////////////////////////// COME FANNO A SPARIRE I PRIMI LATI
    cout<<"Lati: "<<endl;
    for (int i=0;i<int(lati.size());i++)
    {
        Point p1 ;
        p1 = *lati[i].points[0];
        Point p2 ;
        p2 = *lati[i].points[1];

        cout<<to_string(i+1)<<") "<<"punto di inizio: ("<<p1.x<<","<<p1.y<<")       punto di fine: ("<<p2.x<<","<<p2.y <<")"<<endl;

    }


}






/// FUNZIONI





void InputFile(string inputFilePath, vector<Point> &points)
    {
        ifstream file;
        file.open(inputFilePath);
        if (file.fail()){
            cerr<< "file open failed"<< endl;}

        else{
            cout<<"import ok"<<endl;;
        }
        string line;
        //creo una lista di punti
        getline(file,line);
        while (getline(file, line))//while file >> file
        {
         istringstream converter(line); //converter
        int id;
        double x,y;
        converter>>id >> x>>y;
        points.push_back(Point(x,y)); // creo le istanze Punto nella lista points
        }
    file.close();
}



Lati * Point::isconnect(Point*p) // verifico se è collegato con p
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


void Distanzamax(vector<Point> points,Point& p1,Point& p2,Point& p3)
{
 double max_dist = 0.0;
 double max_altezza = 0.0;
 vector<double> distanze;

 // faccio la mergesort per avere algoritmo più efficiente

 for (int i=0;i<int(points.size());i++)  // ciclo sul vettore point che contiene i punti
 {
    for (int j=0;j<int(points.size());j++)  // doppio ciclo
    {
        if (j!=i)
        {
            double distance=Point::dist(points[i],points[j]);
            distanze.push_back(distance);
            if(distance>max_dist)
            {
                max_dist = distance;
                p1 = points[i];
                p2 = points[j];

            }

        }
    }

 }
 double m = ((p2.y - p1.y)/(p2.x - p1.x));
  double q = -m*p1.x+ p1.y;
  for (int i=0;i<int(points.size());i++)
  {
     if (points[i]!= p1 && points[i]!= p2)
     {
         double d = (abs(points[i].y - (m*points[i].x+q))/(sqrt(1+m*m)));
         if (d > max_altezza)
         {
             max_altezza = d;
             p3 = points[i];
         }
     }
  }
 }




bool dentro ( Point a, Point b, Point c)
{
    Point AB=b-a;
    Point AC=c-a;
    double z =AB.x*AC.y-AB.y*AC.x;
    if (z<0)
        return false;
    else
        return true;
}


 bool ispointin(Triangle t, Point &p)
 {
    if(!isInsideCircumcircle(t,p))
        return false;
    else
    {

        // trovo i due punti più vicini a b
        // FUNZIONE MATTE

        Point * a= t.points[0];
        Point * b= t.points[1];
        Point * c= t.points[2];


        if(dentro( *a, * b,  p)&& dentro( *b, * c,  p)&& dentro( *c, * a,  p))
        {
            return true;
        }
        else
            return false;
    }
 }


 bool isInsideCircumcircle(Triangle t, Point& q)
 {
     Point a=*t.points[0];
     Point b=*t.points[1];
     Point c=*t.points[2];

     double ax = a.x - q.x;
     double ay = a.y - q.y;
     double bx = b.x - q.x;
     double by = b.y - q.y;
     double cx = c.x - q.x;
     double cy = c.y - q.y;

     double det = ax * (by * (cx*cx+cy*cy) - cy*(bx*bx+by*by)) - bx * (ay*(cx*cx+cy*cy) - cy*(ax*ax+ay*ay)) + cx * (ay*(bx*bx+by*by) - by*(ax*ax+ay*ay));

     return det > 0.0;

 }


bool siintersecano(Point p1, Point p2, Point p3, Point p4)
{

    int det = (p2.x-p1.x)*(p4.y-p3.y)-(p2.y-p1.y)*(p4.x-p3.x);
    if (det==0)//caso rette parallele
        return false;
    else
    {

        // retta p1p2
        double m1 = ((p2.y - p1.y)/(p2.x - p1.x));
        if (isnan(m1))// retta parallela all'asse x
            m1=10e10;
        double q1 = -m1*p1.x+ p1.y;
        //retta p3p4
        double m2 = ((p4.y - p3.y)/(p4.x - p3.x));
        if (isnan(m2))// retta parallela all'asse x
            m2=10e10;
        double q2 = -m2*p3.x+ p3.y;
        // componente x del punto di intersezione
        double px=(q1-q2)/(m2-m1);
        if (isnan(px))
            px=(p3.x+p4.x)/2;
        // componente y del punto di intersezione
        double py=m1*((q1-q2)/(m2-m1))+q1;
        if (isnan(py))
            py=(p3.y+p4.y)/2;

        Point p(px,py);
        // alpha per la componente x
        double ax=(px-p1.x)/(p2.x-p1.x);
        if (isnan(ax))
            ax=0.00001;
        // alpha per la componente y
        double ay=(py-p1.y)/(p2.y-p1.y);
        if (isnan(ay))
            ay=0.00001;
        // gli alpha devono appartenere a [0,1] per far parte del segmento
        if ((ax>0&&ax<1)&&(ay>0&&ay<1))
            return true;
        else
            return false;

    }

}




void creotriangoli(vector<Triangle> &triangles,int j,Point* p,vector<Lati> &lati)//,vector<Lati> &lati
{
    Triangle old=triangles[j];
    triangles[j]=Triangle(old.points[0],old.points[1],p);
    triangles.push_back(Triangle(old.points[1],old.points[2],p));
    triangles.push_back(Triangle(old.points[2],old.points[0],p));

    lati.push_back(Lati(old.points[0],p,1));
    lati.push_back(Lati(old.points[1],p,1));
    lati.push_back(Lati(old.points[2],p,1));

    // aggiorno le adiacenza dei triangoli
    triangles[j].triangoli_adiacenti.push_back(&triangles[j+1]);
    triangles[j].triangoli_adiacenti.push_back(&triangles[j+2]);
    triangles[j+1].triangoli_adiacenti.push_back(&triangles[j]);
    triangles[j+1].triangoli_adiacenti.push_back(&triangles[j+2]);
    triangles[j+2].triangoli_adiacenti.push_back(&triangles[j+1]);
    triangles[j+2].triangoli_adiacenti.push_back(&triangles[j]);

    //aggoirno la presenza del lato nel triangolo
    for (int i=0;i<3;i++)
    {
        for(int k=0;k<3;k++)
        {
            triangles[j+i].sides[k]->triangoli.push_back(&triangles[j+i]);
        }
    }

}

void Posizione(vector<Point*>& point_da_aggiungere, vector<Triangle> &triangles,vector<Lati*>& lati_di_bordo, vector<Lati>& lati, vector<Point*>& punti_di_bordo)
{
    for (int i=0;i<int(point_da_aggiungere.size());i++)// per ogni punto
    {
        bool in=false;

        for (int j=0;j<int(triangles.size());j++) // per ogni triangolo
        {
            if(ispointin(triangles[j],*point_da_aggiungere[i])) // se il punto sta dentro il triangolo (dereferenzio (*))
            {
                // creo 3 nuovi triangoli
                creotriangoli(triangles,j,point_da_aggiungere[i],lati);//,lati
                in=true;
                ///////////////////////////////////////////////////////////////////////////
                //lati.push_back(*triangles[j].sides[2]);
                //lati.push_back(*triangles[j+1].sides[2]);
                //lati.push_back(*triangles[j+2].sides[2]);

                break; // esco dal ciclo for
            }

        }

        if(!in) // se il punto è esterno a tutti
        {
        vector<Lati *> nuovilati;//vettore con i nuovi lati che non intersecano

           for (int j=0;j<int(punti_di_bordo.size());j++)
           {
                bool intersezione=false;
                Lati l;
                l= Lati(point_da_aggiungere[i],punti_di_bordo[j],1); // creo un nuovo lato

                //Lati* pl;
                //pl=new Lati(point_da_aggiungere[i],punti_di_bordo[j]);
                //cout<<" l puntatore "<<l.points<<endl;
                //cout<<" p1 x "<<l.points[0]->x<<endl;
                //cout<<" p1 y "<<l.points[0]->y<<endl;
                //cout<<" p2 x "<<l.points[1]->x<<endl;
                //cout<<" p2 y "<<l.points[1]->y<<endl;

                //for (int k=0;k<int(lati_di_bordo.size());k++) // verifico che l non intersechi nessun lato esterno
                for (int k=0;k<int(lati.size());k++) // verifico che l non intersechi nessun lato esterno
                {
                    //cout<<" l "<<l.points<<" l2 "<<lati_di_bordo[k]->points<<endl;
                    //cout<<" p1 x "<<lati_di_bordo[k]->points[0]->x<<endl;
                    //cout<<" p1 y "<<lati_di_bordo[k]->points[0]->y<<endl;
                    //cout<<" p2 x "<<lati_di_bordo[k]->points[1]->x<<endl;
                    //cout<<" p2 y "<<lati_di_bordo[k]->points[1]->y<<endl;

                    //if(siintersecano(*l.points[0],*l.points[1],*lati_di_bordo[k]->points[0],*lati_di_bordo[k]->points[1]))
                    if(siintersecano(*l.points[0],*l.points[1],*lati[k].points[0],*lati[k].points[1]))
                    {
                        // se l interseca almeno un lato esterno
                        intersezione=true;
                    }

                }
                //////////////// AGGIUNGE SOLO UN LATO SE IL PUNTO è FUORI

                if(!intersezione) // se non interseca nulla aggiungo il lato
                {
                    lati.push_back(l); // aggiungo il nuovo lato alla lista dei lati esistenti
                    //Lati* pl;
                    //pl=&l;
                    lati_di_bordo.push_back(&l);
                    nuovilati.push_back(&l);
                }



           }


           for(int j=0;j<int(nuovilati.size());j++) // per tutti i nuovi lati creati (validi)
           {
               for (int k=j;k<int(nuovilati.size());k++)
               {
                   if(nuovilati[j]->points[1]->isconnect(nuovilati[k]->points[1])) // trovo il terzo punto del trinagolo che voglio creare
                   {
                       //cout<<" ENTRO NELL IF !!!!!!!!!!! ";
                       Triangle t=Triangle(nuovilati[j]->points[1],nuovilati[k]->points[1],point_da_aggiungere[i]);
                       triangles.push_back(t);
                       //aggiornare triangoli_adiacenti
                       Lati altrolato=Lati(nuovilati[j]->points[1],nuovilati[k]->points[1],1);
                       // aggiungo la adiacenza a quelli esisstenti
                       for(int h=0;h<int(altrolato.triangoli.size());h++) // altrolato.triangoli[0] posso averne solo 2
                       {
                           altrolato.triangoli[h]->triangoli_adiacenti.push_back(&t);
                           t.triangoli_adiacenti.push_back(altrolato.triangoli[h]);
                       }
                       for (int z=0;z<3;z++)
                       {
                           t.sides[z]->triangoli.push_back(&t); // aggiungo alla lista di triangoli di cui l fa parte
                       }

                   }
               }
           }

        punti_di_bordo.push_back(point_da_aggiungere[i]); // aggiungo alla lista di punti di bordo
        }
    ipotesi(point_da_aggiungere[i],lati);
    }

}



void ipotesi(Point* punto,vector<Lati>&lati)
{
    for(int j=0;j<int(punto->edges.size());j++)// per ognilato di cui fa parte
    {
        for (int k=0;k<int(punto->edges[j]->triangoli.size());k++) // per ogni trinagolo con quel lato
        {
            for(int z=0;z<int(punto->edges[j]->triangoli[k]->triangoli_adiacenti.size());z++) // per ogni triangolo adiacente a un triangolo di cui fa parte
            {
                ///////////////////////////////////////// edges 6 non esiste
                if(isInsideCircumcircle(*punto->edges[j]->triangoli[k]->triangoli_adiacenti[z], *punto))
                {
                    Lati oldlato;
                    oldlato=*punto->edges[j];
                    //////////////////////////////////////////////// ERRORE NEL SALVATAGGIO DELLE ADIACENZE!!!!

                    //Flip(*punto->edges[j]->triangoli[k]->triangoli_adiacenti[z],*punto->edges[j]->triangoli[k],oldlato);
                    //remove(lati.begin(), lati.end(), oldlato);

                    //lati.pop_back(); // elimino il lato

                    //Lati nuovo;
                    /// DEVO AGGIUNGERE IL NUOVO LATO
                    //Lati nuovo=*punto->edges[j]->triangoli[k]->triangoli_adiacenti[z]->sides[0];
                    //lati.push_back(nuovo);

                }
            }
        }
    }

}




void Flip(Triangle& ABC, Triangle& CBD, Lati& oldlato) {


   Lati newlato;

   int k = 0;
   for (int i=0; i<3;i++) {                // punti triangolo ABC
       for (int j=0; j<3; j++) {           // punti triangolo CBD
           if (ABC.points[i] == CBD.points[j])
           {
               oldlato.points[k] = ABC.points[i];     //salvo in oldlato il lato adiacente
               k++;

           }
           // se segui il ragionamento del ciclo non dovrebbe sovrascrivere:
           //  A = C, A = B, A = D,
           //  B = C, B = B (SI), B = D,
           //  C = C (SI), C = B, C = D,   -> oldato è points[0] = B, points[1] = C
        }
    }
    //lati.pop_back();

    for (int j=0; j<3; j++){

        if (oldlato.points[0] != ABC.points[j]&&oldlato.points[1] != ABC.points[j])
        {   // prendo il vertice opposto in ABC al lato adiacente
            newlato.points[0] =  ABC.points[j];    // lo salvo in newlato[0]

        };

    }

    for (int j=0; j<3; j++)
    {
        if (oldlato.points[0] != CBD.points[j]&&oldlato.points[1] != CBD.points[j])
        {   // prendo il vertice opposto in CBD al lato adiacente
            newlato.points[1] =  CBD.points[j];    // lo salvo in newlato[1]
        };

    // chiamando i costuttore triangolo a partire dal nuovo punto di newlato crea i triangoli
    }
    // sovrascrivo i triangoli
    ABC = Triangle(newlato.points[0], newlato.points[1], oldlato.points[0] ); //ADB
    CBD = Triangle(newlato.points[0], newlato.points[1], oldlato.points[1] ); //ADC

}

void output (vector<dalaunay::Lati> latiOutput)
{
    int lunghezza = sizeof(latiOutput);
    string a = "";
    cout<<"Lati: "<<endl;
    for (int i=0;i<lunghezza;i++)
    {
        Point p1 ;
        p1 = *latiOutput[i].points[0];
        Point p2 ;
        p2 = *latiOutput[i].points[1];

        cout<<to_string(i+1)<<") "<<"punto di inizio: ("<<p1.x<<","<<p1.y<<")       punto di fine: ("<<p2.x<<","<<p2.y <<")"<<endl;
        //a = to_string(i+1) + ") punto di inizio: (" + to_string(p1.x) + to_string(p1.y)+")       punto di fine: (" + to_string(p2.x)+","+to_string(p2.y)+")";

    }

}













/*

 /// da mettere in ispointin
Point p_giusto1;
Point p_giusto2;

void FunzioneMatte(*p0, *p1, *p2, *p3, &p_giusto1, &p_giusto2)
{

    vector<double> m[3];
    vector<double> q[3];
    vector<double> d[3];
    vector<Point> p[3];
    p.push_back(p1,p2,p3);


    for (int i=0; i<3; i++)
    {
        m[i] = ((p[(i+1)%3].y() - p[i].y())/(p[(i+1)%3].x() - p[i].x()));
        q[i] = -m[i]*p[i].x() + p[i].y();
        d[i] = (abs(p3.y() - (m[i]*p3.x()+q[i]))/(sqrt(1+m[i]*m[i])));
    }

    double minimo=0.0;
    for (int i=0; i<3; i++)
    {
      if (d[i]<minimo)
        minimo=d[i];
    }

    if (minimo==d[2])
    {
        p_giusto1 = &points[0];
        p_giusto2 = &points[2];
    }
    else if (minimo==d[0])
    {
        p_giusto1 = &points[0];
        p_giusto2 = &points[1];
    }
    else
    {
        p_giusto1 = &points[1];
        p_giusto2 = &points[2];
    }
}
*/

/*

double Merge(vector<double> distanze,
           const unsigned int sx,
           const unsigned int cx,
           const unsigned int dx)
{ unsigned int i = sx; unsigned int j = cx+1;
  double massimo = 0.0;
  double massimo_generale = 0.0;
  while ((i <= cx) && (j <= dx))
  {
    if (distanze[i] >= distanze[j])
    {
        massimo = distanze[i]; j = j+1;
    }
    else
    {
      massimo = distanze[j]; i = i+1;
    }
    if (massimo > massimo_generale)
    {
        massimo_generale = massimo;
    }
    }
  return massimo_generale;
   }

double MergeSort(vector<double> distanze,
               const unsigned int sx,
               const unsigned int dx)

{if (sx < dx)
   {
   unsigned int centro = (sx+dx)/2;
   unsigned int meta1 = (sx+centro)/2;
   unsigned int meta2 = (centro+1+dx)/2;
   double a = Merge( distanze, sx, meta1, centro );
   double b = Merge( distanze, centro+1, meta2, dx );
   return (a >b ? a : b);

    }

}


*/



/// DA VERIFICARE


// Verifica l'ipotesi di Delaunay tra due triangoli adiacenti
/*
bool isDelaunayValid( Point& a,  Point& b,  Point& c,  Point& d)
{
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
*/


}
