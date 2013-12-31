#include "assembly.h"
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

/**
 * Assembles the global stiffness matrix given an array of element stiffnesses
 * 
 * @param [mat] kg Pointer to global stiffness matrix
 * @param [MechElem*] elems Array of elements to assemble
 * @param [unsigned int] numElements Number of elements in array
 * @param [unsigned int] dofPerElem Total degrees of freedom per element
 * @return void
 */
void globalStiffness(mat &kg, MechElem **elems, unsigned int numElements,
        unsigned int dofPerElem, int pState)
{
    /*
     * @todo reach a consistency where matrices passed by reference are either \
     * assumed to be zeroed before being passed or are zeroed on pass
     */
    int *gdofs;
    unsigned int elemNum, i, j;
    MechElem *pElem;
    mat *ke = new mat(dofPerElem, dofPerElem);
    for (elemNum = 0; elemNum < numElements; elemNum++)
    {
        pElem = elems[elemNum];
        gdofs = pElem->getGdofs();
        pElem->stiffness(*ke, pState);
        for (i = 0; i < dofPerElem; i++)
        {
            for (j = 0; j < dofPerElem; j++)
            {
                // @todo consider renumbering DOF so that we no longer require \
                // this silly minus one business
                kg(gdofs[i]-1, gdofs[j]-1) += (*ke)(i, j);
            }
        }
    }
    // deallocate dynamic memory
    delete ke;
}