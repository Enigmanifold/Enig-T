function e1s_sph=setup_2Dmesh(k1)
% Construct a densified icosahedron mesh.
    e1s_carte=IcosahedronMesh;
    e1s_carte=SubdivideSphericalMesh(e1s_carte,k1);
    e1s_sph=-ones(2,size(e1s_carte.Points,1));
    e1counter=1;
    for m=1:size(e1s_carte.Points,1)
        if e1s_carte.Points(m,3)>0
            if e1s_carte.Points(m,1)~=0
                sph_coords=[acos(e1s_carte.Points(m,3));atan(e1s_carte.Points(m,2)/e1s_carte.Points(m,1))];
            else
                sph_coords=[acos(e1s_carte.Points(m,3));pi/2];
            end
            if e1s_carte.Points(m,1)<0 && e1s_carte.Points(m,2)>=0
                sph_coords(2)=sph_coords(2)+pi;
            end
            if e1s_carte.Points(m,1)<=0 && e1s_carte.Points(m,2)<0
                sph_coords(2)=sph_coords(2)-pi;
            end
        e1counter=e1counter+1;
        end
        e1s_sph(:,e1counter)=sph_coords;
    end
    e1s_sph=e1s_sph(:,1:e1counter-1);
end