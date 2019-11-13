function cross_mat=cross_matrix(vec)
    vec=vec./norm(vec);
    cross_mat=[0,-vec(3),vec(2);vec(3),0,-vec(1);-vec(2),vec(1),0];
end