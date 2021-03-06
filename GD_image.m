function total_cell=GD_image(p1,p2,ptsnum,maxfunevals,maxiter,densification1,densification2)
    diff_amplifier=1000000;
    e1s_carte=IcosahedronMesh;
    e2s_carte=IcosahedronMesh;
    e1s_carte=SubdivideSphericalMesh(e1s_carte,densification1);
    e2s_carte=SubdivideSphericalMesh(e2s_carte,densification2);
    e1s_sph=[];
    e2s_sph=[];
    counter=1;
    e1counter=1;
    e2counter=1;
    fun=@(x)calc_phi3(x,p1,p2,ptsnum,diff_amplifier);
    options=optimset('MaxFunEvals',maxfunevals,'MaxIter',maxiter,'TolFun',1e-4,'TolX',1e-5);
    total_cell={};
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
            x=[epi1_sph;epi2_sph]';
            min_val_initial=fun(x);
            [xopt,fopt]=fminsearch(fun,x,options);
            total_cell{counter,1}=epi1_sph';
            total_cell{counter,2}=epi2_sph';
            total_cell{counter,3}=xopt(1:2);
            total_cell{counter,4}=xopt(3:4);
            total_cell{counter,5}=min_val_initial;
            total_cell{counter,6}=fopt;
            counter=counter+1;
        end
    end
    total_cell=sortrows(total_cell,6);
end