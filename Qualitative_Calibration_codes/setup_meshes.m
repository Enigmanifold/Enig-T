function epi12_sph=setup_meshes(k1,k2)
% Having densification factors for two icosahedra, create pairs in
% spherical coordinates.
    e1s_carte=IcosahedronMesh;
    e2s_carte=IcosahedronMesh;
    e1s_carte=SubdivideSphericalMesh(e1s_carte,k1);
    e2s_carte=SubdivideSphericalMesh(e2s_carte,k2);
    e1s_sph=-ones(2,size(e1s_carte.Points,1));
    e2s_sph=-ones(2,size(e2s_carte.Points,1));
    e1counter=1;
    e2counter=1;
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
    for m=1:size(e2s_carte.Points,1)
        if e2s_carte.Points(m,3)>0
            if e2s_carte.Points(m,1)~=0
                sph_coords=[acos(e2s_carte.Points(m,3));atan(e2s_carte.Points(m,2)/e2s_carte.Points(m,1))];
            else
                sph_coords=[acos(e2s_carte.Points(m,3));pi/2];
            end
            if e2s_carte.Points(m,1)<0 && e2s_carte.Points(m,2)>=0
                sph_coords(2)=sph_coords(2)+pi;
            end
            if e2s_carte.Points(m,1)<=0 && e2s_carte.Points(m,2)<0
                sph_coords(2)=sph_coords(2)-pi;
            end
            e2counter=e2counter+1;
        end
        e2s_sph(:,e2counter)=sph_coords;
    end
    e2s_sph=e2s_sph(:,1:e2counter-1);
    e1s_sphsize=size(e1s_sph,2);
    e2s_sphsize=size(e2s_sph,2);
    epi12_sph=-ones(e1s_sphsize*e2s_sphsize,4);
    epi12_sph_counter=1;
    for m=1:e1s_sphsize
        for n=1:e2s_sphsize
            epi12_sph(epi12_sph_counter,:)=[e1s_sph(:,m)',e2s_sph(:,n)'];
            epi12_sph_counter=epi12_sph_counter+1;
        end
    end
end