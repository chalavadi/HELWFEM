#include "assembly.h"
#include <armadillo>
#include <iostream>

// Uncomment out the following line to disable assertions
//#define NDEBUG 1
#include <cassert>

using namespace arma;
using namespace std;

/*TODO work with linked lists instead of C-style arrays
 *       this will be easier for memory management */

/**
 * Assembles the global stiffness matrix given an array of elements, mutator
 * 
 * @param [mat] kg Global stiffness matrix (by reference)
 * @param [MechElem*] elems Array of pointers to elements to assemble
 * @param [unsigned int] numElements Number of elements in array
 * @param [unsigned int] dofPerElem Total degrees of freedom per element
 * @return void
 */
void mglobalStiffness(mat &kg, MechElem **elems, unsigned int numElements,
        unsigned int dofPerElem, int pState)
{
    // TODO Can we find a way to relate this to global dof??
    //assert(kg.n_cols == dofPerElem);
    //assert(kg.n_rows == dofPerElem);
    /*
     * TODO reach a consistency where matrices passed by reference are either
     *      assumed to be zeroed before being passed or are zeroed on pass
     */
    kg.zeros();
    int *gdofs;
    unsigned int elemNum, i, j;
    MechElem *pElem;
    mat *ke = new mat(dofPerElem, dofPerElem);
    for (elemNum = 0; elemNum < numElements; elemNum++)
    {
        pElem = elems[elemNum];
        gdofs = pElem->getGdofs();
        pElem->mstiffness(*ke, pState);
        for (i = 0; i < dofPerElem; i++)
        {
            for (j = 0; j < dofPerElem; j++)
            {
                // TODO consider renumbering DOF so that we no longer require
                // this silly minus one business
                kg(gdofs[i]-1, gdofs[j]-1) += (*ke)(i, j);
            }
        }
    }
    // free dynamic memory
    delete ke;
}

/**
 * Assembles the global body force vector given an array of elements, mutator
 * 
 * @param [mat] bg Global body force vector (by reference)
 * @param [MechElem*] elems Array of pointers to elements to assemble
 * @param [unsigned int] numElements Number of elements in array
 * @param [unsigned int] dofPerElem Total degrees of freedom per element
 * @return void
 */
void mglobalBodyForce(vec &bg, MechElem **elems, unsigned int numElements,
        unsigned int dofPerElem)
{
    // TODO Can we write this for global degrees of freedom?
    //assert(bg.n_cols == dofPerElem);
    
    bg.zeros();
    int *gdofs;
    unsigned int elemNum, i;
    MechElem *pElem;
    vec *be = new vec(dofPerElem);
    for (elemNum = 0; elemNum < numElements; elemNum++)
    {
        pElem = elems[elemNum];
        gdofs = pElem->getGdofs();
        pElem->mbodyForce(*be);
        for (i = 0; i < dofPerElem; i++)
            bg(gdofs[i]-1) += (*be)(i);
    }
    // free dynamic memory
    delete be;
}
