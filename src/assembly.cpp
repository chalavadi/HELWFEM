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
void globalStiffness(mat *kg, MechElem **elems, unsigned int numElements,
        unsigned int dofPerElem, int planeState)
{
    cout << "DOF PER ELEM: " << dofPerElem << " num elem: " << numElements << endl;
    unsigned int elemNum, i, j;
    MechElem *pElem;
    mat *ke = new mat(dofPerElem, dofPerElem);
    for (elemNum = 0; elemNum < numElements; elemNum++)
    {
        cout << "element: " << elemNum << endl;
        pElem = elems[elemNum];
        cout << pElem << endl;
        pElem->stiffness(ke, planeState);
        for (i = 0; i < dofPerElem; i++)
            cout << pElem->gdof[i] << endl;
        cout << *ke << endl;
        for (i = 0; i < dofPerElem; i++)
        {
            for (j = 0; j < dofPerElem; j++)
            {
                cout << i << " " << j << endl << *kg << endl;
                cout << pElem->gdof[i] << " " << pElem->gdof[j] << endl;
                (*kg)(pElem->gdof[i], pElem->gdof[j]) += (*ke)(i, j);
            }
        }
    }
    // deallocate dynamic memory
    delete ke;
}