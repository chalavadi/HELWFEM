#include "gtest/include/gtest/gtest.h"
#include "../src/element.h"
#include "../src/assembly.h"
#include <time.h>
#include <iostream>
#include <cstdlib>
//#include "test.h"

#define GTEST_COLOR 1
#define GTEST_OUTPUT 1

using namespace arma;

namespace {

// Define test fixtures
class Q4Test : public ::testing::Test 
{
protected:

    virtual void SetUp() {
        int E = 29000;
        double v = 0.3, h = 1.0;

        b << 0 << endr 
          << -9800 << endr;
        t << 0 << endr
          << 0 << endr;
        gnodes[0] = 1;
        gnodes[1] = 2;
        gnodes[2] = 5;
        gnodes[3] = 4;
        gcoords << 0 << 0 << endr
                << 5 << 0 << endr
                << 5 << 5 << endr
                << 0 << 5 << endr;
        gdofs[0] = 1;
        gdofs[1] = 2;
        gdofs[2] = 3;
        gdofs[3] = 4;
        gdofs[4] = 9;
        gdofs[5] = 10;
        gdofs[6] = 7;
        gdofs[7] = 8;

        elem = Q4(E, v, h, &b, &t, gnodes, &gcoords, gdofs);
        
        gnodes2[0] = 2;
        gnodes2[1] = 3;
        gnodes2[2] = 6;
        gnodes2[3] = 5;
        gcoords2 << 5 << 0 << endr
                << 10 << 0 << endr
                << 10 << 5 << endr
                << 5 << 5 << endr;
        gdofs2[0] = 3;
        gdofs2[1] = 4;
        gdofs2[2] = 5;
        gdofs2[3] = 6;
        gdofs2[4] = 11;
        gdofs2[5] = 12;
        gdofs2[6] = 9;
        gdofs2[7] = 10;
        
        elem2 = Q4(E, v, h, &b, &t, gnodes2, &gcoords2, gdofs2);
    }

    virtual void TearDown() {
    }

    vec b;
    vec t;
    int gnodes[4], gnodes2[4];
    mat gcoords, gcoords2;
    int gdofs[8], gdofs2[8];
    Q4 elem, elem2;
};

/**
 * Test preprocessor shape functions
 */
/* NOT AVAILABLE BECAUSE SHAPE FUNCTIONS ARE DEFINED IN .CPP FILE.
 * CONSIDER MOVING THEM TO HEADER FILE
TEST(Q4Test, ShapeFunctions)
{
}
*/

/**
 * Test the getGcoords method for functionality and polymorphism
 */
TEST_F(Q4Test, Gcoords)
{
    int i, j;
    mat expgcoords(4,2);
    mat *actgcoords = elem.getGcoords();
    
    expgcoords << 0 << 0 << endr
            << 5 << 0 << endr
            << 5 << 5 << endr
            << 0 << 5 << endr;
    
    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            ASSERT_EQ(expgcoords(i,j), (*actgcoords)(i,j));
    
    // Test that polymorphism works the way we want it to
    MechElem *pElem = &elem;
    actgcoords = pElem->getGcoords();
    
    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            ASSERT_EQ(expgcoords(i,j), (*actgcoords)(i,j));
    
    actgcoords = elem2.getGcoords();
    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            ASSERT_EQ(gcoords2(i,j), (*actgcoords)(i,j));
}

/**
 * Test the getGnodes method for functionality and polymorphism
 */
TEST_F(Q4Test, Gnodes)
{
    int i, expgnodes[4], *actgnodes;
    expgnodes[0] = 1;
    expgnodes[1] = 2;
    expgnodes[2] = 5;
    expgnodes[3] = 4;
    
    actgnodes = elem.getGnodes();
    
    for (i = 0; i < 4; i++)
        ASSERT_EQ(expgnodes[i], actgnodes[i]);
    
    // Test inheritance
    MechElem *pElem = &elem;
    actgnodes = pElem->getGnodes();
    
    for (i = 0; i < 4; i++)
        ASSERT_EQ(expgnodes[i], actgnodes[i]);
    
    actgnodes = elem2.getGnodes();
    
    for (i = 0; i < 4; i++)
        ASSERT_EQ(gnodes2[i], actgnodes[i]);
}

/**
 * Test the getGdofs method for functionality and polymorphism
 */
TEST_F(Q4Test, Gdofs)
{
    int i, expgdofs[8], *actgdofs;
    
    expgdofs[0] = 1;
    expgdofs[1] = 2;
    expgdofs[2] = 3;
    expgdofs[3] = 4;
    expgdofs[4] = 9;
    expgdofs[5] = 10;
    expgdofs[6] = 7;
    expgdofs[7] = 8;
    
    actgdofs = elem.getGdofs();
    
    for (i = 0; i < 8; i++)
        ASSERT_EQ(expgdofs[i], actgdofs[i]);
    
    // Test inheritance
    MechElem *pElem = &elem;
    actgdofs = pElem->getGdofs();
    
    for (i = 0; i < 8; i++)
        ASSERT_EQ(expgdofs[i], actgdofs[i]);
    
    actgdofs = elem2.getGdofs();
    
    for (i = 0; i < 8; i++)
        ASSERT_EQ(gdofs2[i], actgdofs[i]);
}

/**
 * Test the stiffness method against values computed in other programs
 */
TEST_F(Q4Test, ElementStiffnessMatrix)
{
    int i, j;
    mat expStiffness(8,8);
    mat actStiffness(8,8);
    expStiffness << 1.4341 << 0.5179 << -0.8764 << -0.0398 << -0.7170 << -0.5179 << 0.1593 << 0.0398 << endr
                  << 0.5179 << 1.4341 << 0.0398 << 0.1593 << -0.5179 << -0.7170 << -0.0398 << -0.8764 << endr
                  << -0.8764 << 0.0398 << 1.4341 << -0.5179 << 0.1593 << -0.0398 << -0.7170 << 0.5179 << endr
                  << -0.0398 << 0.1593 << -0.5179 << 1.4341 << 0.0398 << -0.8764 << 0.5179 << -0.7170 << endr
                  << -0.7170 << -0.5179 << 0.1593 << 0.0398 << 1.4341 << 0.5179 << -0.8764 << -0.0398 << endr
                  << -0.5179 << -0.7170 << -0.0398 << -0.8764 << 0.5179 << 1.4341 << 0.0398 << 0.1593 << endr
                  << 0.1593 << -0.0398 << -0.7170 << 0.5179 << -0.8764 << 0.0398 << 1.4341 << -0.5179 << endr
                  << 0.0398 << -0.8764 << 0.5179 << -0.7170 << -0.0398 << 0.1593 << -0.5179 << 1.4341 << endr;
    expStiffness *= 10000;
    elem.stiffness(actStiffness, PSTRESS);
    
    std::cout << std::endl
              << "Expected and actual stiffness matrices: " << std::endl
              << "----------------------------------------" << std::endl
              << expStiffness << std::endl
              << actStiffness << std::endl
              << "----------------------------------------" << std::endl
                                                            << std::endl;
    
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
            EXPECT_NEAR(expStiffness(i,j), actStiffness(i,j), 
                    abs(expStiffness(i,j)/100));
    
    actStiffness.zeros();
    elem.stiffness(actStiffness, PSTRESS);
    
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
            EXPECT_NEAR(expStiffness(i,j), actStiffness(i,j), 
                    abs(expStiffness(i,j)/100));
    
    expStiffness.zeros();
    expStiffness << 1.4341 <<  0.5179 << -0.8764 << -0.0398 << -0.7170 << -0.5179 <<  0.1593 <<  0.0398 << endr
                 << 0.5179 <<  1.4341 <<  0.0398 <<  0.1593 << -0.5179 << -0.7170 << -0.0398 << -0.8764 << endr
                 << -0.8764 <<  0.0398 <<  1.4341 << -0.5179 <<  0.1593 << -0.0398 << -0.7170 <<  0.5179 << endr
                 << -0.0398 <<  0.1593 << -0.5179 <<  1.4341 <<  0.0398 << -0.8764 <<  0.5179 << -0.7170 << endr
                 << -0.7170 << -0.5179 <<  0.1593 <<  0.0398 <<  1.4341 <<  0.5179 << -0.8764 << -0.0398 << endr
                 << -0.5179 << -0.7170 << -0.0398 << -0.8764 <<  0.5179 <<  1.4341 <<  0.0398 <<  0.1593 << endr
                 << 0.1593 << -0.0398 << -0.7170 <<  0.5179 << -0.8764 <<  0.0398 <<  1.4341 << -0.5179 << endr
                 << 0.0398 << -0.8764 <<  0.5179 << -0.7170 << -0.0398 <<  0.1593 << -0.5179 <<  1.4341 << endr;
    
    expStiffness *= 10000;
    elem2.stiffness(actStiffness, PSTRESS);
    
    std::cout << std::endl
              << "Expected and actual stiffness matrices: " << std::endl
              << "----------------------------------------" << std::endl
              << expStiffness << std::endl
              << actStiffness << std::endl
              << "----------------------------------------" << std::endl
                                                            << std::endl;
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
            EXPECT_NEAR(expStiffness(i,j), actStiffness(i,j), 
                    abs(expStiffness(i,j)/100));
    
    elem2.stiffness(actStiffness, PSTRESS);
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
            EXPECT_NEAR(expStiffness(i,j), actStiffness(i,j), 
                    abs(expStiffness(i,j)/100));
}

TEST_F(Q4Test, GlobalStiffness)
{
    mat ke1(8,8,fill::zeros), ke2(8,8,fill::zeros);
    mat kg(12,12,fill::zeros);
    
    elem.stiffness(ke1, PSTRESS);
    elem2.stiffness(ke2, PSTRESS);
    
    MechElem *pElems[2] = {&elem, &elem2};
    
    globalStiffness(kg, pElems, 2, 8, PSTRESS);
    
    std::cout << std::endl << "Global stiffness matrix" << std::endl
         << "----------------------------------------" << std::endl
         << kg << std::endl;
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_DOUBLE_EQ(ke1(i,j), kg(i,j));
    
    for (int i = 2; i < 4; i++)
        for (int j = 2; j < 4; j++)
            EXPECT_DOUBLE_EQ(ke1(i,j)+ke2(i-2,j-2), kg(i,j));
    
    for (int i = 4; i < 6; i++)
        for (int j = 4; j < 6; j++)
            EXPECT_DOUBLE_EQ(ke2(i-2,j-2), kg(i,j));
    
    for (int i = 6; i < 8; i++)
        for (int j = 6; j < 8; j++)
            EXPECT_DOUBLE_EQ(ke1(i,j), kg(i,j));
    
    for (int i = 8; i < 10; i++)
        for (int j = 8; j < 10; j++)
            EXPECT_DOUBLE_EQ(ke1(i-4,j-4)+ke2(i-2,j-2),kg(i,j));
    
    for (int i = 10; i < 12; i++)
        for (int j = 10; j < 12; j++)
            EXPECT_DOUBLE_EQ(ke2(i-6,j-6), kg(i,j));
}

TEST_F(Q4Test, ElementBodyForce)
{
    vec expBodyForce(8);
    expBodyForce << 0
                 << -61250
                 << 0
                 << -61250
                 << 0
                 << -61250
                 << 0
                 << -61250;
    vec actBodyForce(8, fill::zeros);
    elem.bodyForce(actBodyForce);
    for (int i = 0; i < 8; i++)
        EXPECT_DOUBLE_EQ(expBodyForce(i),actBodyForce(i));
    elem2.bodyForce(actBodyForce);
    for (int i = 0; i < 8; i++)
        EXPECT_DOUBLE_EQ(expBodyForce(i),actBodyForce(i));
    
    std::cout << std::endl << "Element body force" << std::endl
         << "----------------------------------------" << std::endl
         << actBodyForce << std::endl;
}

TEST_F(Q4Test, GlobalBodyForce)
{
    vec expBodyForce(12);
    expBodyForce << 0
                 << -61250
                 << 0
                 << -61250*2
                 << 0
                 << -61250
                 << 0
                 << -61250
                 << 0
                 << -61250*2
                 << 0
                 << -61250;
    vec actBodyForce(12, fill::zeros), be1(8, fill::zeros), be2(8, fill::zeros);
    elem.bodyForce(be1);
    elem2.bodyForce(be2);
    MechElem *pElems[2] = {&elem, &elem2};
    globalBodyForce(actBodyForce, pElems, 2, 8);
    for (int i = 0; i < 12; i++)
        EXPECT_DOUBLE_EQ(expBodyForce(i),actBodyForce(i));
    std::cout << std::endl << "Global body force" << std::endl
         << "----------------------------------------" << std::endl
         << actBodyForce << std::endl;
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}