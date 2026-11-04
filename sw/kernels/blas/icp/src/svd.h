#include <stdlib.h>
#include <math.h>

void jacobi_SVD_3x3(float matrix[9], float U[9]){
    int maxiter = 20;
    float tolerance = 10e-6f;

    for (int iter = 0; iter < maxiter; iter++) {
        // initialize U as the identity-matrix
        // float U = {1, 0, 0,
        //           0, 1, 0,
        //           0, 0, 1};
                   
        // find largest off-diagonal value and save it's coordinates
        int i_coord = 0;
        int j_coord = 1;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (!(i == j)) {
                    if (fabs(matrix[i * 3 + j]) > fabs(matrix[i_coord * 3 + j_coord])) {
                        i_coord = i;
                        j_coord = j;   
                    }
                }
            }
        }

        if (fabs(matrix[i_coord * 3 + j_coord]) < tolerance) {
            break;
        }


        // calculate tan(2a)
        float Sii = matrix[3 * i_coord + i_coord];
        float Sjj = matrix[3 * j_coord + j_coord];
        float Sij = matrix[3 * i_coord + j_coord];
        float cotan2a = (Sjj - Sii) / (2 * Sij);
        // calculate tan(a)
        float tana;
        if (cotan2a >= 0) {
            tana = 1.0f / (cotan2a + sqrt(1.0f + cotan2a * cotan2a));
        } else {
            tana = -1.0f / (-cotan2a + sqrt(1.0f + cotan2a * cotan2a));
        }
        // calculate cos(a) and sin(a)
        float cosa = 1.0f / (sqrt(1.0f + tana * tana));
        float sina = tana / (sqrt(1.0f + tana * tana));

        // update matrix_copy
        matrix[3 * i_coord + i_coord] = cosa * cosa * Sii - 2.0f * sina * cosa * Sij + sina * sina * Sjj;
        matrix[3 * j_coord + j_coord] = sina * sina * Sii + 2.0f * sina * cosa * Sij + cosa * cosa * Sjj;
        matrix[3 * i_coord + j_coord] = 0.0f;
        matrix[3 * j_coord + i_coord] = 0.0f;
        for (int k = 0; k < 3; k++) {
            if (k != i_coord && k != j_coord) {
                float mjk = matrix[3 * j_coord + k];
                float mik = matrix[3 * i_coord + k];
                matrix[3 * i_coord + k] = cosa * mik - sina * mjk;
                matrix[3 * k + i_coord] = matrix[3 * i_coord + k];

                matrix[3 * j_coord + k] = sina * mik + cosa * mjk;
                matrix[3 * k + j_coord] = matrix[3 * j_coord + k];
            }
        }

        // populate U
        for (int k = 0; k < 3; k++) {
            float Uki = U[3 * k + i_coord];
            float Ukj = U[3 * k + j_coord];
            U[3 * k + i_coord] = cosa * Uki - sina * Ukj;
            U[3 * k + j_coord] = sina * Uki + cosa * Ukj;
        }
    }
}