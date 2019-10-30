gte1s_carte=IcosahedronMesh;
gte2s_carte=IcosahedronMesh;
% k1=1;
% k2=1;
% gte1s_carte=SubdivideSphericalMesh(gte1s_carte,k1);
% gte2s_carte=SubdivideSphericalMesh(gte2s_carte,k2);
gte1s_sph=[];
gte2s_sph=[];
gte1counter=1;
gte2counter=1;
for m=1:size(gte1s_carte.Points,1)
    if gte1s_carte.Points(m,3)>0
        if gte1s_carte.Points(m,1)~=0
            sph_coords=[acos(gte1s_carte.Points(m,3));atan(gte1s_carte.Points(m,2)/gte1s_carte.Points(m,1))];
        else
            sph_coords=[acos(gte1s_carte.Points(m,3));pi/2];
        end
        if gte1s_carte.Points(m,1)<0 && gte1s_carte.Points(m,2)>=0
            sph_coords(2)=sph_coords(2)+pi;
        end
        if gte1s_carte.Points(m,1)<=0 && gte1s_carte.Points(m,2)<0
            sph_coords(2)=sph_coords(2)-pi;
        end
    gte1counter=gte1counter+1;
    end
    gte1s_sph(:,gte1counter)=sph_coords;
end
for m=1:size(gte2s_carte.Points,1)
    if gte2s_carte.Points(m,3)>0
        if gte2s_carte.Points(m,1)~=0
            sph_coords=[acos(gte2s_carte.Points(m,3));atan(gte2s_carte.Points(m,2)/gte2s_carte.Points(m,1))];
        else
            sph_coords=[acos(gte2s_carte.Points(m,3));pi/2];
        end
        if gte2s_carte.Points(m,1)<0 && gte2s_carte.Points(m,2)>=0
            sph_coords(2)=sph_coords(2)+pi;
        end
        if gte2s_carte.Points(m,1)<=0 && gte2s_carte.Points(m,2)<0
            sph_coords(2)=sph_coords(2)-pi;
        end
        gte2counter=gte2counter+1;
    end
    gte2s_sph(:,gte2counter)=sph_coords;
end
gte1s_sphsize=size(gte1s_sph,2);
gte2s_sphsize=size(gte2s_sph,2);
gtepis12_sph=-ones(gte1s_sphsize*gte2s_sphsize,4);
gtepis12_sphcounter=1;
for m=1:gte1s_sphsize
    for n=1:gte2s_sphsize
        gtepis12_sph(gtepis12_sphcounter,:)=[gte1s_sph(:,m)',gte2s_sph(:,n)'];
        gtepis12_sphcounter=gtepis12_sphcounter+1;
    end
end