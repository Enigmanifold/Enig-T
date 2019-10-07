run('generate_RT.m');
ptsnum=30;
e1s_carte=IcosahedronMesh;
e2s_carte=IcosahedronMesh;
k1=3;
k2=3;
e1s_carte=SubdivideSphericalMesh(e1s_carte,k1);
e2s_carte=SubdivideSphericalMesh(e2s_carte,k2);
e1s_sph=[];
e2s_sph=[];
for m=1:length(size(e1s_carte.Points,1))
    if e1s_carte.Points(m,3)>0
        sph_coords=[acos(e1s_carte.Points(m,3));acos(e1s_carte.Points(m,2)/e1s_carte.Points(m,1))];
        if e1s_carte.Points(m,2)<0
            sph_coords(2)=-sph_coords(2);
        end
    end
    e1s_sph(:,end+1)=sph_coords;
end
for m=1:length(size(e2s_carte.Points,1))
    if e2s_carte.Points(m,3)>0
        sph_coords=[acos(e2s_carte.Points(m,3));acos(e2s_carte.Points(m,2)/e2s_carte.Points(m,1))];
        if e2s_carte.Points(m,2)<0
            sph_coords(2)=-sph_coords(2);
        end
    end
    e2s_sph(:,end+1)=sph_coords;
end
for m=1:size(e1s_sph,2)
    for n=1:size(e2s_sph,2)
        epi1_sph=e1s_sph(:,m);
        epi2_sph=e2s_sph(:,n);
        epi1_plane=[tan(epi1_sph(1))*cos(epi1_sph(2));tan(epi1_sph(1))*sin(epi1_sph(2)),1];
        epi2_plane=[tan(epi2_sph(1))*cos(epi2_sph(2));tan(epi2_sph(1))*sin(epi2_sph(2)),1];
        [phi1,phi2,phi1_alt,phi2_alt]=calc_phi(p1,p2,epi1_plane,epi2_plane,ptsnum);
        