#include "assembly.h"
#include <armadillo>

using namespace arma;

/**
 * Assembles the global stiffness matrix given an array of element stiffnesses
 * 
 * @param [mat] kg Pointer to global stiffness matrix
 * @param [MechElem*] elems Array of elements to assemble
 * @param [unsigned int] numElements Number of elements in array
 * @param [unsigned int] dofPerElem Total degrees of freedom per element
 * @return void
 */
void globalStiffness(mat *kg, MechElem *elems, unsigned int numElements,
        unsigned int dofPerElem, int planeState)
{
    unsigned int elemNum, i, j;
    MechElem *pElem;
    mat ke = new mat(dofPerElem, dofPerElem);
    for (pElem = elems, elemNum = 0; elemNum < numElements; pElem++, elemNum++)
    {
        pElem->stiffness(&ke, planeState);
        for (i = 0; i < dofPerElem; i++)
        {
            for (j = 0; j < dofPerElem; j++)
            {
                *kg(pElem->gdof[i], pElem->gdof[j]) += ke(i, j);
            }
        }
    }
    // deallocate dynamic memory
    delete ke;
}