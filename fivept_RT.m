function [fiveptR,fiveptT]=fivept_RT(p1_pixel,p2_pixel,ptsnum,K)
% Algorithm parameters
params = [];
params.RANSACIterations = 5000;                         % 5000 iterations
params.reprojectionErrorThreshold = 0.001;                  % 2 pixels
params.bidirectionalMatchConsistencyThreshold = 0.001;      % 2 pixels
% params.correspondenceDensificationFactor = 4; 
p1_pixel=p1_pixel(:,1:ptsnum);
p2_pixel=p2_pixel(:,1:ptsnum);
f1_match = p1_pixel(1:2,:);
f2_match = p2_pixel(1:2,:);
f2_match_rows = p2_pixel'; % Each row represents one pix position.
f1_match_cols = p1_pixel; % Each column represents one pix position.
%Find E
best_inliers = -1;
best_E = zeros(3,3);
best_inds = [];
invK=inv(K);
for iter = 1:params.RANSACIterations
    candi_inds = randsample(ptsnum,5)';
    pixels1 = f1_match(:,candi_inds);
    pixels2 = f2_match(:,candi_inds);
    Es = ComputeEssentialMatrix(pixels1,pixels2,K);
    for i = 1:size(Es,3)
        E = Es{i};
        errs = zeros(1,ptsnum);
        F = invK'*E*invK;
        abc = F*f1_match_cols;
        ab = abc(1:2,:);
        norm_const = vecnorm(ab,2,1);
        for n = 1:ptsnum
            err = abs(f2_match_rows(n,:)*abc(:,n));
            errs(1,n) = err;
        end
        errs = errs./norm_const;
        good_inds = find(errs < params.reprojectionErrorThreshold);
        inliers = length(good_inds);
        if inliers > best_inliers
            best_inliers = inliers;
            best_E = E;
            best_inds = good_inds;
        end
    end
end
% Drawing epipolar lines.
best_F = invK'*best_E*invK;
best_inliers = length(best_inds);
f1_matches = f1_match(:,best_inds);
f2_matches = f2_match(:,best_inds);
% xy1 = round(f1_matches(:,2800));
% xyz1 = [xy1;1];
% % abc_ = best_F*xyz1;
% abc_ = xyz1'*best_F;
% a = abc_(1);
% b = abc_(2);
% c = abc_(3);
% xline = 1:size(image2,2);
% slope = -a/b;
% intersect = -c/b;
% yline = (-a/b.*xline)-(c/b).*ones(size(xline));
% figure
% imshow(image1)
% hold on
% plot(xy1(1),xy1(2),'*r');
% hold off
% figure
% imshow(image2)
% hold on
% plot(xline,yline,'r');
% hold off
% Calculate R and T.
f2_match_rows = [f2_matches',ones(best_inliers,1)]; % Each row represents one pix position.
f1_match_cols = [f1_matches;ones(1,best_inliers)]; % Each column represents one pix position.
EE_trans = best_E*best_E';
TT_trans = 0.5.*trace(EE_trans)*eye(3)-EE_trans;
T = TT_trans(1,:)'./sqrt(TT_trans(1,1));
Tx = [0,-T(3),T(2);T(3),0,-T(1);-T(2),T(1),0];
R = 1/(T'*T)*(det(best_E).*inv(best_E))'-1/(T'*T).*Tx*best_E;
minusT = -T;
rotR = 1/(minusT'*minusT)*(det(best_E).*inv(best_E))'-1/(minusT'*minusT).*(-Tx)*best_E;
votes = zeros(1,4);
gamma1 = zeros(3,ptsnum);
gamma2 = zeros(3,ptsnum);
for m = 1:best_inliers
    g1 = K\f1_match_cols(:,m);
    g2 = K\f2_match_rows(m,:)';
    gamma1(:,m) = g1;
    gamma2(:,m) = g2;
end
for n = 1:best_inliers
    rg1 = -R*gamma1(:,n);
    rg2 = -rotR*gamma1(:,n);
    b1 = T(1:3);
    b2 = minusT(1:3);
    A1 = [gamma2(1,n),rg1(1);gamma2(2,n),rg1(2);gamma2(3,n),rg1(3)];
    A2 = [gamma2(1,n),rg2(1);gamma2(2,n),rg2(2);gamma2(3,n),rg2(3)];
    rhos1 = A1\b1;    
    rhos2 = A1\b2;
    rhos3 = A2\b1;
    rhos4 = A2\b2;
    rhos = [rhos1,rhos2,rhos3,rhos4];
    for v = 1:4
        if rhos(1,v) > 0 && rhos(2,v) > 0
            votes(v) = votes(v)+1;
        end
    end
end
switch max(votes)
    case votes(1)
        fiveptR = R;
        fiveptT = T;
    case votes(2)
        fiveptR = R;
        fiveptT = minusT;
    case votes(3)
        fiveptR = rotR;
        fiveptT = T;
    case votes(4)
        fiveptR = rotR;
        fiveptT = minusT;
end