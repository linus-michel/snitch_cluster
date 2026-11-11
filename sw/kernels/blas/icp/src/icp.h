void icp_single_core(float *ref_pc, float *src_pc, int N)
{
    // ref_pc: reference point cloud
    // src_pc: source point cloud
    // N: amount of points in the clouds

    // TODO: maybe make that this is passed in the function
    int maxiter = 5;

    for (int iterations = 0; iterations < 5; iterations++)
    {
        // calculate centroids
        float centroid_ref_pc[3] = {0.0f, 0.0f, 0.0f};
        float centroid_src_pc[3] = {0.0f, 0.0f, 0.0f};

        // temp
        // [0, 0, 0, 0]
        // [1, 1, 1, 1]
        // [2, 2, 2, 2]

        // [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]

        // need to drop homogeneous coordinates for the following computations
        // sum all points
        for (int i = 0; i < N; i++)
        {
            centroid_ref_pc[0] += ref_pc[0 * N + i];
            centroid_ref_pc[1] += ref_pc[1 * N + i];
            centroid_ref_pc[2] += ref_pc[2 * N + i];

            centroid_src_pc[0] += src_pc[0 * N + i];
            centroid_src_pc[1] += src_pc[1 * N + i];
            centroid_src_pc[2] += src_pc[2 * N + i];
        }

        // divide by N to get center of the point clouds
        for (int i = 0; i < 3; i++)
        {
            centroid_ref_pc[i] /= N;
            centroid_src_pc[i] /= N;
        }

        // center the point clouds
        float ref_pc_centered[192];
        float src_pc_centered[192];
        for (int i = 0; i < 192; i++)
        {
            ref_pc_centered[i] = ref_pc[i];
            src_pc_centered[i] = src_pc[i];
        }

        for (int i = 0; i < N; i++)
        {
            ref_pc_centered[0 * N + i] -= centroid_ref_pc[0];
            ref_pc_centered[1 * N + i] -= centroid_ref_pc[1];
            ref_pc_centered[2 * N + i] -= centroid_ref_pc[2];

            src_pc_centered[0 * N + i] -= centroid_src_pc[0];
            src_pc_centered[1 * N + i] -= centroid_src_pc[1];
            src_pc_centered[2 * N + i] -= centroid_src_pc[2];
        }

        // temp
        // [3]  x   [0, 1, 2]
        // [4]
        // [5]

        // [3x0, 3x1, 3x2]
        // [4x0, 4x1, 4x2]
        // [5x0, 5x1, 5x2]

        // calculate H = sum(src*ref^T)   -   fixed bug that it didn't calculate H = sum(src*ref^T) but H = sum(ref*src^T)
        float H[9] = {0.0f};
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                H[3 * j + 0] += ref_pc_centered[0 * N + i] * src_pc_centered[j * N + i];
                H[3 * j + 1] += ref_pc_centered[1 * N + i] * src_pc_centered[j * N + i];
                H[3 * j + 2] += ref_pc_centered[2 * N + i] * src_pc_centered[j * N + i];
            }
        }

        // temp
        // [0, 3, 6]        [0, 1, 2]       [0*0 + 3*3 + 6*6, 0*1 + 3*4 + 6*7, 0*2 + 3*5 + 6*8]         [45, 54, 63]             H =     [0, 1, 2]
        // [1, 4, 7]    x   [3, 4, 5]   =   [1*0 + 4*3 + 7*6, 1*1 + 4*4 + 7*7, 1*2 + 4*5 + 7*8]     =   [54, 66, 78]                     [3, 4, 5]
        // [2, 5, 8]        [6, 7, 8]       [2*0 + 5*3 + 8*6, 2*1 + 5*4 + 8*7, 2*2 + 5*5 + 8*8]         [63, 78, 93]                     [6, 7, 8]

        
        // calculate S = HTH    -   works
        float S[9];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                S[j * 3 + i] = (H[0 + i] * H[0 + j] + H[3 + i] * H[3 + j] + H[6 + i] * H[6 + j]);
            }
        }
        
        // initialize V to the identity matrix
        float V[9] = {1, 0, 0,
                      0, 1, 0,
                      0, 0, 1};

        // perform SVD to get V, since V of H are the eigenvectors of HTH
        jacobi_SVD_3x3(S, V);

        // order the eigenvalues by size in descending order:
        // first find the right order
        int a = 0;
        int b = 1;
        int c = 2;
        int t = 0;

        // ensures a > b
        if (S[3 * a + a] < S[3 * b + b]) {
            t = a;
            a = b;
            b = t;
        }
        // ensures a > b > c
        if (S[3 * a + a] < S[3 * c + c]) {
            t = c;
            c = b;
            b = a;
            a = t;
        }
        // ensures b > c
        else if (S[3 * b + b] < S[3 * c + c]) {
            t = c;
            c = b;
            b = t;
        }
        


        // order S and V by the previously found order
        // 2 0 0    a = 1, b = 0, c = 2     =>  change columns 0 and 1 of V
        // 0 3 0
        // 0 0 1

        int order[3] = {a, b, c};
        float V_ordered[9];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                V_ordered[3 * j + i] = V[j * 3 + order[i]];
            }
        }

        // initialize U to the identity matrix
        float U[9] = {1, 0, 0,
                      0, 1, 0,
                      0, 0, 1};
        
        // calculate U by using Ui = 1/Si * HVi
        for (int j = 0; j < 3; j++) {
            float sigma = sqrt(S[order[j] * 3 + order[j]]);
            for (int i = 0; i < 3; i++) {
                // TODO: divide by 0 issues if sigma is small, problem for later
                U[3 * i + j] = (1 / sigma) * (H[3 * i + 0] * V_ordered[3 * 0 + j] + H[3 * i + 1] * V_ordered[3 * 1 + j] + H[3 * i + 2] * V_ordered[3 * 2 + j]);
            }
        }


        // calculate X=VU^T
        float X[9];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                X[3 * i + j] = V_ordered[3 * i + 0] * U[3 * j + 0] + V_ordered[3 * i + 1] * U[3 * j + 1] + V_ordered[3 * i + 2] * U[3 * j + 2];
            }
        }

        // check if det(X) = 1
        float X_determinant = X[0] * (X[4] * X[8] - X[5] * X[7]) -
                              X[1] * (X[3] * X[8] - X[5] * X[6]) +
                              X[2] * (X[3] * X[7] - X[4] * X[6]);

        float rotation_matrix[9];
        if (X_determinant > 0.0f)
        {
            for (int i = 0; i < 9; i++)
            {
                rotation_matrix[i] = X[i];
            }
        }
        else
        {
            V_ordered[2] = V_ordered[2] * -1;
            V_ordered[5] = V_ordered[5] * -1;
            V_ordered[8] = V_ordered[8] * -1;

           for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                X[3 * i + j] = V_ordered[3 * i + 0] * U[3 * j + 0] + V_ordered[3 * i + 1] * U[3 * j + 1] + V_ordered[3 * i + 2] * U[3 * j + 2];
            }
        }
        }

        // calculate translation vector
        // first apply the rotation matrix to the centroid of the source pointcloud
        float rotated_src_centroid[3];
        for (int i = 0; i < 3; i++)
        {
            rotated_src_centroid[i] = rotation_matrix[3 * i + 0] * centroid_src_pc[0] +
                                      rotation_matrix[3 * i + 1] * centroid_src_pc[1] +
                                      rotation_matrix[3 * i + 2] * centroid_src_pc[2];
        }

        // then subtract rotated_src_pc from centroid_ref_pc
        float translation_vector[3];
        translation_vector[0] = centroid_ref_pc[0] - rotated_src_centroid[0];
        translation_vector[1] = centroid_ref_pc[1] - rotated_src_centroid[1];
        translation_vector[2] = centroid_ref_pc[2] - rotated_src_centroid[2];

        // construct transformation matrix
        // R = rotation matrix, T = translation vector
        //
        // R R R T
        // R R R T
        // R R R T
        // 0 0 0 1
        float transformation_matrix[16] = {
            rotation_matrix[0], rotation_matrix[1], rotation_matrix[2], translation_vector[0],
            rotation_matrix[3], rotation_matrix[4], rotation_matrix[5], translation_vector[1],
            rotation_matrix[6], rotation_matrix[7], rotation_matrix[8], translation_vector[2],
            0.0f,               0.0f,               0.0f,               1.0f
        };


        // temp
        //  trans_matrix        qi
        //  m0  m1  m2  m3      v0          m0v0 + m1v1 + m2v2 + m3v3
        //  m4  m5  m6  m7      v1          m4v0 + m5v1 + m6v2 + m7v3
        //  m8  m9  m10 m11     v2          ...
        //  m12 m13 m14 m15     v3          ...

        // apply transformation
        float transformed_pc[256];
        for (int i = 0; i < N; i++)
        {
            float src_0ni = src_pc[0 * N + i];
            float src_1ni = src_pc[1 * N + i];
            float src_2ni = src_pc[2 * N + i];
            float src_3ni = src_pc[3 * N + i];

            transformed_pc[0 * N + i] = transformation_matrix[0] * src_0ni + 
                                        transformation_matrix[1] * src_1ni + 
                                        transformation_matrix[2] * src_2ni + 
                                        transformation_matrix[3] * 1 /*src_3ni*/;
            transformed_pc[1 * N + i] = transformation_matrix[4] * src_0ni + 
                                        transformation_matrix[5] * src_1ni + 
                                        transformation_matrix[6] * src_2ni + 
                                        transformation_matrix[7] * 1 /*src_3ni*/;
            transformed_pc[2 * N + i] = transformation_matrix[8] * src_0ni + 
                                        transformation_matrix[9] * src_1ni + 
                                        transformation_matrix[10] * src_2ni + 
                                        transformation_matrix[11] * 1 /*src_3ni*/;
            /*
            transformed_pc[3 * N + i] = transformation_matrix[3] * src_0ni + 
                                        transformation_matrix[7] * src_1ni + 
                                        transformation_matrix[11] * src_2ni + 
                                        transformation_matrix[15] * src_3ni;
            */
           // TODO: fix w coordinate calculation
           transformed_pc[3 * N + i] = src_3ni;
        }


        // scale the x,y,z such that w matches the w from the input
        /* for (int i = 0; i < N; i++) {
            float w = src_pc[3 * N + i];
            float w_trans = transformed_pc[3 * N + i];

            if (fabs(w_trans) > 1e-3) {
                float scaling = w / w_trans;
                transformed_pc[0 * N + i] = transformed_pc[0 * N + i] / scaling;
                transformed_pc[1 * N + i] = transformed_pc[1 * N + i] * scaling;
                transformed_pc[2 * N + i] = transformed_pc[2 * N + i] * scaling;
                transformed_pc[3 * N + i] = w;
            }
        }  */


        // put the transformed_pc back into src_pc
        for (int i = 0; i < 256; i++) {
            src_pc[i] = transformed_pc[i];
        }

        // convergence check
        // TODO: implement convergence check
        if (iterations >= maxiter) {
            break;
        }
    }
}