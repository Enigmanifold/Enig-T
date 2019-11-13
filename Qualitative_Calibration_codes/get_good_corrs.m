function [p1_pixel,p2_pixel,p1,p2,good_corrs]=get_good_corrs(im1,im2,R_gt,T_gt,K1,K2)
% Have two images, find correspondences that give small reprojection
% errors.
    params = [];
    params.ratioTestRatio = 0.75;
    params.RANSACIterations = 1000;
    params.RANSACInlierDistanceThreshold = 0.5;
    [good_pairs,f1,f2] = find_good_pairs(im1,im2,params);
    pts_size=size(good_pairs,2);
    p1_pixel = [f1(1:2,good_pairs(1,:));ones(1,pts_size)];
    p2_pixel = [f2(1:2,good_pairs(2,:));ones(1,pts_size)];
    p1=K1\p1_pixel;
    p2=K2\p2_pixel;
    good_corrs=zeros(1,pts_size);
    good_corrs_counter=1;
    max_reproj_err=0.2/max(K1(1),K2(1));
    for m=1:pts_size
        re=sum(compute_reprojection_error_all(R_gt,T_gt,p1(:,m),p2(:,m),1))/size(p1,2);
        if re<max_reproj_err
            good_corrs(good_corrs_counter)=m;
            good_corrs_counter=good_corrs_counter+1;
        end
    end
    good_corrs=good_corrs(1:good_corrs_counter-1);
end