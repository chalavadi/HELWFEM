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
    
  vec be1 = zeros<colvec>(8);
  vec be2 = zeros<colvec>(8);
  
  int gnodes[] = {1, 2, 5, 4};
  mat gcoords(4,2);
  gcoords << 0 << 0 << endr
          << 5 << 0 << endr
          << 5 << 5 << endr
          << 0 << 5 << endr;
  int gdofs[] = {1, 2, 3, 4, 9, 10, 7, 8};
  
  Q4 *elem = new Q4(E, v, h, &b, &t, gnodes, &gcoords, gdofs);
  
  int gnodes2[] = {2, 3, 6, 5};
  mat gcoords2(4,2);
  gcoords2 << 0 << 0 << endr
          << 5 << 0 << endr
          << 5 << 5 << endr
          << 0 << 5 << endr;
  int gdof2[] = {3, 4, 5, 6, 11, 12, 9, 10};
  
  Q4 *elem2 = new Q4(E, v, h, &b, &t, gnodes2, &gcoords2, gdof2);
  
  mat ke(8,8);
  elem->stiffness(ke, PSTRESS);
  
  mat ke2(8,8);
  elem2->stiffness(ke2, PSTRESS);
  
  MechElem *pelems[2];
  pelems[0] = elem;
  pelems[1] = elem2;
  
  mat kg = zeros<mat>(12, 12);
  
  int numElem = 2;
  globalStiffness(kg, pelems, numElem, Q4__DOF_PER_NODE*Q4__NUM_NODES, PSTRESS);
  
  cout << "Element Stiffness 1" << endl
       << "-------------" << endl
       << ke << endl
       << "Element Stiffness 2" << endl
       << "-------------" << endl
       << ke2 << endl
       << "Global Stiffness" << endl
       << "-------------" << endl
       << kg << endl;
  
  elem->bodyForce(be1);
  elem2->bodyForce(be2);
  
  vec bg(12);
  globalBodyForce(bg, pelems, numElem, Q4__DOF_PER_NODE*Q4__NUM_NODES);
  
  cout << endl << "Body Force 1" << endl
       << "-------------" << endl
       << be1 << endl
       << "-------------" << endl
       << "Body Force 2" << endl
       << be2 << endl
       << "Global Body Force" << endl
       << bg << endl;

  delete elem;
  delete elem2;
  
  return 0;
}
