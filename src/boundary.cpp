#include <list>

//#define NDEBUG 1
#include <cassert>
#include "boundary.h"

/**
 * Impose boundary conditions for a model, mutator
 *
 * @param kg Global stiffness matrix passed by reference
 * @param fg Global force vector passed by reference
 * @param bounds List of boundary conditions
 */
void mimposeBoundaryConds(mat &kg, vec &fg, const std::list<BC*> &bounds)
{
    assert(kg.n_rows == kg.n_cols);
    assert(kg.n_rows == fg.n_rows);
    
    int i;
    int gdof = kg.n_rows;
    for (auto bc : bounds)
    {
        for (i = 0; i < gdof; i++)
        {
            if (i == bc->dof)
            {
                kg(i, bc->dof) = 1;
                continue;
            }
            /* modify force vector with known products */
            fg(i) -= kg(i, bc->dof) * bc->val - kg(bc->dof, i) * bc->val;
            /* zero off-diagonals */
            kg(i, bc->dof) = 0;
            kg(bc->dof, i) = 0;
        }
    }
    
    for (auto bc: bounds)
    /* modify force vector with known displacements
     * this is done after original loop to ensure none of the EBCs are
     *      over written */
        fg(bc->dof) = bc->val;
}
