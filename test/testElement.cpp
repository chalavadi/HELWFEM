#include "gtest/include/gtest/gtest.h"
#include "../src/element.h"
#include <time.h>
#include <iostream>
#include <cstdlib>
//#include "test.h"

#define GTEST_COLOR 1
#define GTEST_OUTPUT 1

using namespace arma;

namespace {

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
        gnodes[2] = 4;
        gnodes[3] = 5;
        gcoords << 0 << 0 << endr
                << 5 << 0 << endr
                << 5 << 5 << endr
                << 0 << 5 << endr;
        gdof[0] = 1;
        gdof[1] = 2;
        gdof[2] = 3;
        gdof[3] = 4;
        gdof[4] = 9;
        gdof[5] = 10;
        gdof[6] = 7;
        gdof[7] = 8;

        elem = Q4(E, v, h, &b, &t, gnodes, &gcoords, gdof);
    }

    virtual void TearDown() {
    }

    vec b;
    vec t;
    int gnodes[4];
    mat gcoords;
    int gdof[8];
    Q4 elem;
};

TEST_F(Q4Test, ElementStiffnessMatrix) {
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
    elem.stiffness(&actStiffness, PSTRESS);
    
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
    actStiffness = elem.stiffness(PSTRESS);
    
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
            EXPECT_NEAR(expStiffness(i,j), actStiffness(i,j), 
                    abs(expStiffness(i,j)/100));
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}