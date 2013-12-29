#ifndef TEST_H
#define	TEST_H

#include <cstdlib>

#define EXPECT_AEQ_PER(exp, act, per) EXPECT_LE(abs(exp-act), abs(exp*per/100))

/*
#include <cstdlib>
#include <cassert>
#include <time.h>

#define EPS_SMALL 0.001
#define EPS_LARGE 10

#define CHECK(expression, msg) \
if (!(expression)) \
{ \
    passed = 0; \
    std::cout << msg << endl; \
}

#define CHECK_EQ(exp, act) CHECK(exp == act)
#define CHECK_ALEQ_EPS(exp, act, eps) CHECK(abs(exp - act) < eps)
#define CHECK_ALEQ_PER(exp, act, per) CHECK(abs(exp - act) < exp * per/100)
#define CHECK_LT(exp, act) CHECK(exp < act)
#define CHECK_GT(exp, act) CHECK(exp > act)
#define CHECK_LE(exp, act) CHECK(exp <= act)
#define CHECK_GE(exp, act) CHECK(exp >= act)
#define CHECK_NE(exp, act) CHECK(exp != act)

#define INIT_SUITE \
time_t suiteStart = time(NULL); \
unsigned int testsRan = 0; \
int passed;

#define TEST(name) \
std::cout << name \
passed = 1;

#define END_TEST \
if (!(passed)) \
{ \
    std::cout << "FAILED" << endl;
*/

#endif	/* TEST_H */