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
  int gnodes[] = {1, 2, 5, 4};
  mat gcoords(4,2);
  gcoords << 0 << 0 << endr
          << 5 << 0 << endr
          << 5 << 5 << endr
          << 0 << 5 << endr;
  int gdofs[] = {1, 2, 3, 4, 9, 10, 7, 8};
  
  Q4 *elem = new Q4(E, v, h, &b, &t, gnodes, &gcoords, gdofs);
  
  /*
  int gnodes2[] = {2, 3, 6, 5};
  mat gcoords2(4,2);
  gcoords2 << 0 << 0 << endr
          << 5 << 0 << endr
          << 5 << 5 << endr
          << 0 << 5 << endr;
  int gdof2[] = {3, 4, 5, 6, 11, 12, 9, 10};
  
  Q4 *elem2 = new Q4(E, v, h, &b, &t, gnodes2, &gcoords2, gdof2);
  
  mat ke(8,8);
  elem->stiffness(&ke, PSTRESS);
  
  mat ke2(8,8);
  elem2->stiffness(&ke2, PSTRESS);
  
  MechElem *pelems[2];
  pelems[0] = elem;
  pelems[1] = elem2;
  
  mat::fixed<12,12> kg;
  kg.zeros();
  
  int numElem = 2;
  globalStiffness(&kg, pelems, numElem, Q4__DOF_PER_NODE*Q4__NUM_NODES, PSTRESS);
  
  cout << "Element 1" << endl
       << "-------------" << endl
       << ke << endl
       << "Element 2" << endl
       << "-------------" << endl
       << ke2 << endl
       << "Global" << endl
       << "-------------" << endl
       << kg << endl;
  */
  
  mat *pgcoords;
  pgcoords = elem->getGcoords();
  
  cout << *pgcoords << endl;
  
  int *new_nodes;
  new_nodes = elem->getGnodes();
  for (int i = 0; i < 4; i++)
      cout << new_nodes[i] << endl;
  
  int *new_gdof;
  new_gdof = elem->getGdofs();
  for (int i = 0; i < 8; i++)
      cout << new_gdof[i] << endl;

  delete elem;
  
  return 0;
}