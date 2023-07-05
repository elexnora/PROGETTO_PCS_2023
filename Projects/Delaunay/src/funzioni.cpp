#include "funzioni.hpp"



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
    void LeggoFileInput(const string &inputFilePath, vector<Point>& punti)
        {
            ifstream file; // apro file in modalità lettura
            file.open(inputFilePath);
            if (file.fail())
            {
                cerr<< "file open failed"<< endl;
            }

            else
            {
                cout<<"import ok"<<endl;
            }
            string line;
            getline(file,line);  // leggo ogni linea del file
            int id = 0;
            while (getline(file,line))//while file >> file
            {
             istringstream converter(line); //converter
             double x,y;
             converter>>id>>x>>y;
             punti.push_back(Point(x,y,id)); // inserisco il nuovo punto letto nella lista punti
             id++;// creo le istanze Punto nella lista points
            }
        file.close();
    }

    void CreoFileOutput(const string &outputFilePath,const  vector<Lati> &myList)
    {
        ofstream outputFile(outputFilePath); // Apriamo un nuovo file in modalità di scrittura

        if (outputFile.is_open())
        {
            outputFile<<"X1 Y1 X2 Y2"<<endl;
            for (const auto& element : myList) // per ogni elemento nella lista
            {
                Point p1 ;
                p1 = element.points[0];
                Point p2 ;
                p2 = element.points[1];

                outputFile << p1.x  << " " << p1.y << " " << p2.x << " " << p2.y << "\n"; // Scrivo ogni elemento della lista separato da uno spazio nel file, mettendo separate le componenti delle x e delle y dei punt finali e iniziali del lato
            }

            outputFile.close(); // Chiudiamo il file dopo aver terminato la scrittura
            cout << "Lista scritta nel file con successo!" << std::endl;
        }
        else
        {
            cout << "Impossibile aprire il file." << endl;
        }

    }

}
