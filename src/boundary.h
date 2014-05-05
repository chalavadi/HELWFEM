#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <armadillo>
#include <list>

using namespace std;
using namespace arma;

/**
 * Boundary condition
 *
 * @param dof Degree of freedom
 * @param val Known value
 */
class BC
{
public:
    int dof;
    double val;
    BC();
    BC(const int dof, const double val) : dof(dof), val(val) {};
};


void mimposeBoundaryConds(mat&, vec&, const list<BC*>&);

#endif /* BOUNDARY_H */
