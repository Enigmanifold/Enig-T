function [R,T,p1,p2]=episph_phi_to_RT(epi12_sph,phi)
    [epi1_plane(1),epi1_plane(2)]=imsphere2implane(epi12_sph(1),epi12_sph(2));
    [epi2_plane(1),epi2_plane(2)]=imsphere2implane(epi12_sph(3),epi12_sph(4));
    epi1_plane=[epi1_plane,1]';
    epi2_plane=[epi2_plane,1]';
    flag = 0;
    while (flag == 0)
        T_scale=randsample(2,1)*2-3;
        [R1,R2,T1,T2]=epipoles_phi_to_RT(epi1_plane,epi2_plane,phi,T_scale);
        rnd1=rand;
        rnd2=rand;
        if rnd1 < 0.5
            R=R1;
        else
            R=R2;
        end
        if rnd2 < 0.5
            T=T1;
        else
            T=T2;
        end
        %R = eul2rotmatrix(eulRot,'ZYX');
        % intrinsic matrix K is pre-defined & Also the image size is pre-defined.
        % K =   [718.8560         0  607.1928;
        %          0  718.8560  185.2157;
        %          0         0    1.0000;];
        K = [250 0 101;
             0 250 101;
             0 0 1];
        %K = eye(3);
        imSize = [101*2, 101*2];

        % Build Essential Matrix tx * R
        tx = [0 -T(3) T(2);
              T(3)  0 -T(1);
              -T(2) T(1) 0];

        EGT = tx * R;

        %%Generate 3D Points 
        DX = (rand(1,20000) - 0.5) * 5;
        DY = (rand(1,20000) - 0.5) * 5;
        DZ = (rand(1,20000) + 1.5) * 2;
        DPoints = [DX; 
                   DY; 
                   DZ];

        %%Projected on images
        for i = 1:size(DPoints,2)
            p1(:,i) = DPoints(:,i) ./ DPoints(3,i);
            p1Image(:,i) = K * p1(:,i);
        end

        for i = 1:size(DPoints,2)
            P3D2(:,i) = R * (DPoints(:,i)) + T;
            p2(:,i) = P3D2(:,i) ./ P3D2(3,i);
            p2Image(:,i) = K * p2(:,i);
        end

        index1X = p1Image(1,:) < imSize(1) & p1Image(1,:) > 1;
        index1Y = p1Image(2,:) < imSize(1) & p1Image(2,:) > 1;
        index2X = p2Image(1,:) < imSize(1) & p2Image(1,:) > 1;
        index2Y = p2Image(2,:) < imSize(1) & p2Image(2,:) > 1;
        index3 = P3D2(3,:) > 0;

        pointIndex = index1X & index1Y & index2X & index2Y & index3;

        DPointsOnImage = DPoints(:,pointIndex);
        p1ImageOnImage = p1Image(:,pointIndex);
        p2ImageOnImage = p2Image(:,pointIndex);

        if size(DPointsOnImage,2) > 40
            flag = 1;
        end
    end
end