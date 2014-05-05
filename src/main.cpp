#include <iostream>
#include <armadillo>
#include <list>

#include "element.h"
#include "assembly.h"
#include "boundary.h"

using namespace std;
using namespace arma;

//TODO create material.h
//TODO parse input file to initialize model data
//TODO parse command line args
//TODO it may make more sense to restructure data, each element containing
//  pointers to node objects?

int main(int argc, char** argv)
{
  cout << "Program launched." << endl;
  int i;

  /* initialize model data */
  const double h = 0.1;
  const long long int E = 200e9;
  const double v = 0.3;
  const int rho = 7800;
  const double g = 9.8;
  
  /* body and traction forces */
  vec b(Q4__DOF_PER_NODE);
  b     << 0            << endr 
        << -(rho * g)   << endr;
  vec t(Q4__DOF_PER_NODE);
  t     << 0            << endr
        << 0            << endr;
  
  /* boundary conditions */
  list<BC*> bounds;
  bounds.push_front(new BC(0, 0));
  bounds.push_front(new BC(1, 0));
  bounds.push_front(new BC(3, 0));
  bounds.push_front(new BC(5, 0));
  
  /* mesh */
  cout << "Creating mesh..." << endl;
  
  const int numElem = 2;
  const int numNodes = 6;
  
  int gnodes[] = {1, 2, 5, 4};
  mat gcoords(Q4__NUM_NODES, Q4__DOF_PER_NODE);
  gcoords   << 0.    << 0.    << endr
            << 10.   << 2.    << endr
            << 5.    << 8.    << endr
            << 2.    << 6.    << endr;
  int gdofs[] = {1, 2, 3, 4, 9, 10, 7, 8};
  Q4 *elem = new Q4(E, v, h, &b, &t, gnodes, &gcoords, gdofs);
  
  int gnodes2[] = {2, 3, 6, 5};
  mat gcoords2(Q4__NUM_NODES, Q4__DOF_PER_NODE);
  gcoords2  << 10.   << 2.    << endr
            << 20.   << 0.    << endr
            << 17.   << 6.    << endr
            << 5.    << 8.    << endr;
  int gdof2[] = {3, 4, 5, 6, 11, 12, 9, 10};
  Q4 *elem2 = new Q4(E, v, h, &b, &t, gnodes2, &gcoords2, gdof2);
  
  /* calculate element stiffnesses and assemble */
  cout << "Analyzing discretized system..." << endl;
  
  MechElem *pelems[numElem] = {elem, elem2};
  
  mat kg = zeros<mat>(Q4__DOF_PER_NODE*numNodes, Q4__DOF_PER_NODE*numNodes);
  mglobalStiffness(kg, pelems, numElem, Q4__DOF_PER_NODE*Q4__NUM_NODES, PSTRESS);
  
  cout << endl << "Global Stiffness Matrix" << endl
       << kg << endl;
  
  /* calculate element body forces and assemble */
  vec bg(Q4__DOF_PER_NODE*numNodes);
  mglobalBodyForce(bg, pelems, numElem, Q4__DOF_PER_NODE*Q4__NUM_NODES);
  
  cout << endl << "Global Body Force" << endl
       << bg << endl;
  
  /* assemble global force vector */
  vec fg = zeros<vec>(Q4__DOF_PER_NODE*numNodes);
  fg(4) = 10.e3;
  fg(7) = -10.e3;
  fg(9) = -10.e3;
  for (i = 0; i < Q4__DOF_PER_NODE*numNodes; i++)
    fg(i) += bg(i);
    
  cout << endl << "Global Force" << endl
       << fg << endl;
    
  /* impose boundary condtions */
  cout << endl << "Imposing boundary conditions..." << endl;
  mimposeBoundaryConds(kg, fg, bounds);
  
  /* solve for displacements */
  cout << "Solving for displacements..." << endl;
  vec ug = zeros<vec>(Q4__DOF_PER_NODE*numNodes);
  if (!solve(ug, kg, fg))
  {
    cout << endl << "Error: solution not found" << endl;
    return 1;
  }
  
  /* display answer */
  cout << endl << "================================" << endl;
  cout << endl << "Modified Global Stiffness Matrix (by imposing EBCs)" << endl
       << kg << endl;
  cout << endl << "Modified Global Force (by imposing EBCs)" << endl
       << fg << endl;
  cout << endl << "Global Displacements" << endl
       << ug << endl;

  /* clean up dynamic memory */
  cout << endl << "Cleaning up allocated memory..." << endl;
  for (BC *pbc : bounds) delete pbc;
  delete elem;
  delete elem2;
  
  return 0;
}
