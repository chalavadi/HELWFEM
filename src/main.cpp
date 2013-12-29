#include <iostream>
#include <armadillo>
#include "element.h"
#include "assembly.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
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
  
  mat ke(8,8);
  elem->stiffness(&ke, PSTRESS);
  cout << ke << endl;
  cout << trans(ke) << endl;
  
  delete elem;
  
  return 0;
}