#ifndef __TEST_EMPTY_COPY_H
#define __TEST_EMPTY_COPY_H

#include <gtest/gtest.h>

#include "empty_class.hpp"


using namespace testing;
using namespace dalaunay;

TEST(TestDelaunay, TestOrdinaPunti)
{

    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);

    Triangle Triangolo(&p1,&p2,&p3);
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(Triangolo, TriangoloOrdinato);
}

TEST(TestDelaunay,TestIsPointOutside)
{
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point pext(5.0,6.0);
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(ispointin(TriangoloOrdinato,pext), false);
}

TEST(TestDelaunay,TestIsPointOutside2)
{
    Point p1(0.915736,0.138624);
    Point p2(0.0971318,0.917194) ;
    Point p3(0.0357117,0.254282);
    Point pext(0.814724,0.276025) ;
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(ispointin(TriangoloOrdinato,pext), false);
}

TEST(TestDelaunay,TestIsPointInside)
{
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point pint(0.25,0.25);
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(ispointin(TriangoloOrdinato,pint), true);
}

TEST(TestDelaunay,TestIsPointBetween)
{
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point pbetween(0.6,0.6);
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(ispointin(TriangoloOrdinato,pbetween), false);
}

TEST(TestDelaunay,TestMaxDistance)
{

    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point p4(10.0, 10.0);
    Point p5(2.0,0.0);
    Point p6(0.0,2.0);
    Point p7(1.5,0.0);
    Point p8(0.0,1.2);
    Point p9(5.0,20.0);
    Point p10(0.6,0.6);
    Point &a=p5;
    Point &b=p2;
    Point &c=p3;
    vector<Point> punti;
    punti.push_back(p1);
    punti.push_back(p2);
    punti.push_back(p3);
    punti.push_back(p4);
    punti.push_back(p5);
    punti.push_back(p6);
    punti.push_back(p7);
    punti.push_back(p8);
    punti.push_back(p9);
    punti.push_back(p10);
    Distanzamax( punti, a, b, c);
    EXPECT_EQ(a+b+c,p1+p4+p9);
}

TEST(TestDelaunay,TestIntersectionSegments)
{
    Point p1(0.0,0.0);
    Point p2(2.0,0.0);
    Point p3(1.0,-1.0);
    Point p4(1.0,1.0);


    EXPECT_EQ(siintersecano(p1,p2,p3,p4),true);
}
TEST(TestDelaunay,TestIntersectionSegments2)
{
    Point p1(0.0,0.0);
    Point p2(2.0,0.0);
    Point p3(3.0,-1.0);
    Point p4(3.0,1.0);

    EXPECT_EQ(siintersecano(p1,p2,p3,p4),false);
}

TEST(TestDelaunay,TestIntersectionSegments3)
{
    Point p1(0.0,0.0);
    Point p2(2.0,1.0);
    Point p3(1.0,-2.0);
    Point p4(3.0,2.0);

    EXPECT_EQ(siintersecano(p1,p2,p3,p4),false);
}

TEST(TestDelaunay, TestFlip)
{
    Point p1(0.0,0.0);
    Point p2(1.0,0.0);
    Point p3(0.0,1.0);
    Triangle ABC (&p1,&p2,&p3);
    Point p4(2.0, 1.0);
    Triangle CBD (&p3,&p2,&p4);
    Lati lato(&p1,&p2);
    //Lati &l=lato;
    Flip(ABC, CBD,lato);
    Triangle ADB (&p1,&p4,&p2);
    Triangle ADC (&p1,&p4,&p3);
    EXPECT_EQ(ABC,ADB);
    EXPECT_EQ(CBD,ADC);
    EXPECT_EQ(lato,Lati(&p2,&p3));


}

TEST(TestDelaunay,TestCreoTriangoli){
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point pInt(0.25,0.25);

    int j=0;

    Triangle Triangolo(&p1,&p3,&p2);

    vector<Triangle> ListaTriangoli;
    vector<Triangle> ListaTriangoliOutput;
    vector<Lati> lati;

    ListaTriangoli.push_back(Triangolo);

    Triangle Triangolo1(&p1,&p3,&pInt);
    Triangle Triangolo2(&pInt,&p3,&p2);
    Triangle Triangolo3(&p1,&pInt,&p2);

    ListaTriangoliOutput.push_back(Triangolo1);
    ListaTriangoliOutput.push_back(Triangolo2);
    ListaTriangoliOutput.push_back(Triangolo3);
    creotriangoli(ListaTriangoli,j,&pInt,lati);//,lati

    EXPECT_EQ(ListaTriangoli,ListaTriangoliOutput);
}

TEST(TestDelaunay, TestInput)
{
    vector <Point> points;
    string InputFilePath = "C:/Users/matteo/OneDrive/Desktop/PCS2023_Exercises/Projects/Delaunay/Dataset/Test2.csv";
    vector<Point> punti{
    Point(0.814723686393179, 0.276025076998578),
    Point( 0.905791937075619, 0.679702676853675),
    Point( 0.126986816293506, 0.655098003973841),
    Point( 0.913375856139019, 0.162611735194631),
    Point( 0.63235924622541, 0.118997681558377),
    Point( 0.0975404049994095, 0.498364051982143),
    Point( 0.278498218867048, 0.959743958516081),
    Point( 0.546881519204984, 0.340385726666133),
    Point( 0.957506835434298, 0.585267750979777),
    Point( 0.964888535199277, 0.223811939491137),
    Point( 0.157613081677548, 0.751267059305653),
    Point( 0.970592781760616 ,0.255095115459269),
    Point( 0.957166948242946, 0.505957051665142),
    Point( 0.485375648722841, 0.699076722656686),
    Point( 0.8002804688888, 0.890903252535798),
    Point( 0.141886338627215, 0.959291425205444),
    Point( 0.421761282626275, 0.547215529963803),
    Point( 0.915735525189067, 0.138624442828679),
    Point( 0.792207329559554, 0.149294005559057),
    Point( 0.959492426392903, 0.257508254123736),
    Point( 0.655740699156587, 0.840717255983663),
    Point( 0.0357116785741896, 0.254282178971531),
    Point( 0.849129305868777, 0.814284826068816),
    Point( 0.933993247757551, 0.243524968724989),
    Point( 0.678735154857773, 0.929263623187228),
    Point( 0.757740130578333, 0.349983765984809),
    Point( 0.743132468124916, 0.196595250431208),
    Point( 0.392227019534168, 0.251083857976031),
    Point( 0.655477890177557, 0.616044676146639),
    Point( 0.171186687811562, 0.473288848902729),
    Point( 0.706046088019609, 0.351659507062997),
    Point( 0.0318328463774207, 0.830828627896291),
    Point( 0.27692298496089, 0.585264091152724),
    Point( 0.0461713906311539, 0.54972360829114),
    Point( 0.0971317812358475, 0.91719366382981),
    Point( 0.823457828327293, 0.285839018820374),
    Point( 0.694828622975817, 0.757200229110721),
    Point( 0.317099480060861, 0.753729094278495),
    Point( 0.950222048838355, 0.380445846975357),
    Point( 0.0344460805029088, 0.567821640725221),
    Point( 0.438744359656398, 0.0758542895630636),
    Point( 0.381558457093008, 0.0539501186666071),
    Point( 0.765516788149002, 0.530797553008973),
    Point( 0.795199901137063, 0.779167230102011),
    Point( 0.186872604554379, 0.934010684229183),
    Point( 0.489764395788231, 0.12990620847373),
    Point( 0.445586200710899, 0.568823660872193),
    Point( 0.646313010111265, 0.469390641058206),
    Point( 0.709364830858073, 0.0119020695012414),
    Point(0.754686681982361, 0.337122644398882)};

    InputFile(InputFilePath,points);
    EXPECT_EQ(points,punti);
}

TEST(TestDelaunay, TestInput2)
{
    vector <Point> points;
    string InputFilePath = "C:/Users/matteo/OneDrive/Desktop/PCS2023_Exercises/Projects/Delaunay/Dataset/prova.csv";
    vector<Point> punti{
    Point(0.0000000000000000e+00, 0.0000000000000000e+00),
    Point( 1.5000000000000000e+00 ,0.0000000000000000e+00),
    Point( 1.0000000000000000e+00, 1.0000000000000000e+00),
    Point( 0.1000000000000000e+00 ,1.1000000000000000e+00),
    Point( 4.0000000000000000e-01, 5.0000000000000000e-01),
    Point( 1.1000000000000000e+00 ,3.0000000000000000e-01),
    Point( 2.5000000000000000e-01, 0.7000000000000000e+00),
    Point(1.3000000000000000e-00 ,1.5000000000000000e+00),
    Point(9.0000000000000000e-01, 1.10000000000000000e+00),
    Point(4.0000000000000000e-01 ,0.3000000000000000e+00)
    };

    InputFile(InputFilePath,points);
    EXPECT_EQ(points,punti);

}
/*
TEST(TestDelaunay,TestOutput)
{
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point p4(1.0,1.0);
    Point p5(2.0,2.0);
    Point p6(3.4,3.1);
    vector<Lati> lativerifica{ Lati(&p1, &p2), Lati(&p3,&p4), Lati(&p5,&p6)};
    //string b = to_string(1) + ") punto di inizio: (" + to_string(p1.x) + to_string(p1.y)+")       punto di fine: (" + to_string(p2.x)+","+to_string(p2.y)+")";
    output (lativerifica);
    //EXPECT_EQ(b, a);


}
*/


/// PROVA

TEST(TestDelaunay,TestMappa)
{
    Point p1(0.0,0.0);
    Point p2(2.0,1.0);
    Point p3(1.0,-2.0);
    Point p4(3.0,2.0);
    map<int, Point> p;
        p[1] = p1;

        pair<int, Point> elem1 = {2, p2};
        pair<int, Point> elem2 = {3, p3};

        p.insert(elem1);
        p.insert(elem2);
    cout<<p[1].x;
    for(int i=0; i<int(p.size());i++)
    {}
    std::map<int, char> map = { {1, 'A'}, {2, 'B'}, {3, 'C'}, {4, 'D'} };
        int key = 4;

        map.erase(key);

        for (const auto &entry: map) {
            std::cout << "{" << entry.first << ", " << entry.second << "}" << std::endl;
        }
    //EXPECT_EQ(siintersecano(p1,p2,p3,p4),false);
}










#endif // __TEST_EMPTY_COPY_H
