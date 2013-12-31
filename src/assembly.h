#ifndef ASSEMBLY_H
#define	ASSEMBLY_H

#include <armadillo>
#include "element.h"

using namespace arma;

void globalStiffness(mat&, MechElem**, unsigned int, unsigned int, int);
void globalBodyForce(vec&, MechElem**, unsigned int, unsigned int);

#endif	/* ASSEMBLY_H */

