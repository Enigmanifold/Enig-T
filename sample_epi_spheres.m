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
total_cell={};
counter=1;
e1counter=1;
e2counter=1;
e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
if e1(1)<0
    e1Phi=e1Phi+pi;
end
e2Phi=atan(e2(2)/e2(1));
if e2(1)<0
    e2Phi=e2Phi+pi;
end
e1_sph=[e1Theta;e1Phi];
e2_sph=[e2Theta;e2Phi];
e1s_sph(:,1)=e1_sph;
e2s_sph(:,1)=e2_sph;
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
for m=1:size(e1s_sph,2)
    for n=1:size(e2s_sph,2)
        epi1_sph=e1s_sph(:,m);
        epi2_sph=e2s_sph(:,n);
        epi1_angdist=acos(e1sinTheta*sin(epi1_sph(1))+e1cosTheta*cos(epi1_sph(1))*cos(epi1_sph(2)-epi1_sph(2)));
        epi2_angdist=acos(e2sinTheta*sin(epi2_sph(1))+e2cosTheta*cos(epi2_sph(1))*cos(epi2_sph(2)-epi2_sph(2)));
        epi1_plane=[tan(epi1_sph(1))*cos(epi1_sph(2));tan(epi1_sph(1))*sin(epi1_sph(2));1];
        epi2_plane=[tan(epi2_sph(1))*cos(epi2_sph(2));tan(epi2_sph(1))*sin(epi2_sph(2));1];
        [phi1,phi2,phi1_alt,phi2_alt]=calc_phi(p1(:,1:ptsnum),p2(:,1:ptsnum),epi1_plane,epi2_plane,ptsnum);
        phis=[phi1;phi2;phi1_alt;phi2_alt];
        phis_c0=-ones(4,ptsnum);
        for c=1:4               
            for d=1:ptsnum
                [~,phi2_c0_ind]=min([abs(phis(c,d)),abs(phis(c,d)-pi),abs(phis(c,d)-2*pi)]);
                switch phi2_c0_ind
                    case 1
                        phis_c0(c,d)=phis(c,d);
                    case 2
                        phis_c0(c,d)=phis(c,d)-pi;
                    case 3
                        phis_c0(c,d)=phis(c,d)-2*pi;
                end
            end
        end
        phi1_c0std=std(phis_c0(1,:));
        phi2_c0std=std(phis_c0(2,:));
        phi3_c0std=std(phis_c0(3,:));
        phi4_c0std=std(phis_c0(4,:));
        [minc0std,minc0stdind]=min([phi1_c0std,phi2_c0std,phi3_c0std,phi4_c0std]);
        chosen_phis=phis(minc0stdind,:);
        total_cell{counter,1}=epi1_sph';
        total_cell{counter,2}=epi2_sph';
        total_cell{counter,3}=epi1_angdist;
        total_cell{counter,4}=epi2_angdist;
        total_cell{counter,5}=epi1_plane';
        total_cell{counter,6}=epi2_plane';
        total_cell{counter,7}=chosen_phis;
        total_cell{counter,8}=minc0std;
        counter=counter+1;
    end
end