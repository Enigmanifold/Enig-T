function Es = findE(f1,f2,match,K)
    pairs_len = size(match,2);
    f1_pos = f1(1:2,:);
    f2_pos = f2(1:2,:);
%     for iter = 1:params.RANSACIterations
    candi_inds = randsample(pairs_len,4)';
    pixels1 = f1_pos(:,candi_inds);
    pixels2 = f2_pos(:,candi_inds);
    Es = ComputeEssentialMatrix(pixels1,pixels2,K);
end