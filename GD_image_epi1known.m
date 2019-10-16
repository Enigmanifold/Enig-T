function [GD_x,total_cell]=GD_image_epi1known(epi1_sph,p1,p2,ptsnum,maxfunevals,maxiter,densification2)
    diff_amplifier=1000000;
    e2s_carte=IcosahedronMesh;
    e2s_carte=SubdivideSphericalMesh(e2s_carte,densification2);
    e2s_sph=[];
    counter=1;
    e2counter=1;
    fun=@(x)calc_phi3_epi1known(x,epi1_sph,p1,p2,ptsnum,diff_amplifier);
    options=optimset('MaxFunEvals',maxfunevals,'MaxIter',maxiter,'TolFun',1e-4,'TolX',1e-5);
    total_cell=[];
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
    for n=1:size(e2s_sph,2)
        epi2_sph=e2s_sph(:,n);
        x=epi2_sph';
        [xopt,fopt]=fminsearch(fun,x,options);
        total_cell(counter,1:2)=xopt;
        total_cell(counter,3)=fopt;
        counter=counter+1;
    end
    total_cell=sortrows(total_cell,3);
    GD_x=sort(total_cell(1:5,1:2));
    GD_x=GD_x(2,:)+GD_x(3,:)+GD_x(4,:);
    GD_x=GD_x./3;
end