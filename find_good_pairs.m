function [good_pairs,f1,f2,d1,d2] = find_good_pairs(im1,im2,params)
    [f1,d1] = vl_sift(im1);
    [f2,d2] = vl_sift(im2);
    dists = zeros(size(d1,2),size(d2,2));
    % Construct Euclidean distance matrix for all feature pairs.
    for m = 1:size(d1,2)
        for n = 1:size(d2,2)
            displacement = single(d1(:,m)-d2(:,n));
            dists(m,n) = norm(displacement);
        end
    end
    % Ratio test. Obtain indices of matched pairs.
    ind1 = zeros(1,min(size(d1,2),size(d2,2)));
    ind2 = zeros(1,min(size(d1,2),size(d2,2)));
    pair_num = 1;
    s_dists = sort(dists,2);
    for n = 1:size(dists,1)
        if s_dists(n,1) < params.ratioTestRatio*s_dists(n,2)
            ind1(pair_num) = n;
            [~,index2] = min(dists(n,:));
            ind2(pair_num) = index2;
            pair_num = pair_num+1;
        end
    end
    pair_num = pair_num - 1;
    ind1 = ind1(1,1:(pair_num));
    ind2 = ind2(1,1:(pair_num));
    good_pairs = vertcat(ind1,ind2);
end