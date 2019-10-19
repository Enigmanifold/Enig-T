function ang_dist=calculate_ang_dist(vec_sph1,vec_sph2)
    ang_dist=acos(cos(vec_sph1(1))*cos(vec_sph2(1))+sin(vec_sph1(1))*sin(vec_sph2(1))*cos(vec_sph1(2)-vec_sph2(2)));
end