#ifndef T_H
#define T_H

#endif // T_H



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

    EXPECT_EQ(TriangoloOrdinato.ispointin(pext), false);
}

TEST(TestDelaunay,TestIsPointInside)
{
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point pint(0.25,0.25);
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(TriangoloOrdinato.ispointin(pint), true);
}

TEST(TestDelaunay,TestIsPointBetween)
{
    Point p1(0.0,0.0);
    Point p2(0.0,1.0);
    Point p3(1.0,0.0);
    Point pbetween(0.6,0.6);
    Triangle TriangoloOrdinato(&p1,&p3,&p2);

    EXPECT_EQ(TriangoloOrdinato.ispointin(pbetween), false);
}
#endif // __TEST_EMPTY_H
