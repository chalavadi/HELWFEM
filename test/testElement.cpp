#include "gtest/include/gtest/gtest.h"
#include "../src/element.h"
#include <time.h>
#include <iostream>
#include <cstdlib>
//#include "test.h"

#define GTEST_COLOR 1
#define GTEST_OUTPUT 1

namespace {
/*
class Q4Test : public ::testing::Test 
{
protected:
    Q4Test()
    {
        
    }

    virtual ~Q4Test()
    {
        // pass
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp() {
        _startTime = time(NULL);
        
        int E = 29000;
        double v = 0.3, h = 1.0;

        vec b(2);
        b << 0 << endr 
          << -9800 << endr;
        vec t(2);
        t << 0 << endr
          << 0 << endr;
        int gnodes[] = {1, 2, 4, 5};
        mat gcoords(4,2);
        gcoords << 0 << 0 << endr
                << 5 << 0 << endr
                << 5 << 5 << endr
                << 0 << 5 << endr;
        int gdof[] = {1, 2, 3, 4, 9, 10, 7, 8};

        mat *pgcoords = &gcoords;

        elem = new Q4(E, v, h, &b, &t, gnodes, pgcoords, gdof);
    }

    virtual void TearDown() {
        const time_t _endTime = time(NULL);
        
        delete elem;
    }

    Q4 *elem;
    time_t _startTime;
}; */

TEST(Foo, Foobar)
{
    ASSERT_TRUE(1);
}

TEST(Q4Test, ElementStiffnessMatrix)
{
    int E = 29000;
    double v = 0.3, h = 1.0;

    vec b(2);
    b << 0 << endr 
      << -9800 << endr;
    vec t(2);
    t << 0 << endr
      << 0 << endr;
    int gnodes[] = {1, 2, 4, 5};
    mat gcoords(4,2);
    gcoords << 0 << 0 << endr
            << 5 << 0 << endr
            << 5 << 5 << endr
            << 0 << 5 << endr;
    int gdof[] = {1, 2, 3, 4, 9, 10, 7, 8};

    mat *pgcoords = &gcoords;

    Q4 *elem = new Q4(E, v, h, &b, &t, gnodes, pgcoords, gdof);
    
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
    elem->stiffness(&actStiffness, PSTRESS);
    
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
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}