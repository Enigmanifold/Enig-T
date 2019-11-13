run_num=1000;
ptsnum=100;
errs=-ones(run_num,ptsnum);
fiveptRs=cell(run_num,1);
fiveptTs=cell(run_num,1);
for m=1:run_num
    [R_gt,T_gt,e1_gt,e2_gt,phi_gt,K,p1OnImage,p2OnImage,p1_pixel,p2_pixel]=gen_RT;
    p1_pixel=p1_pixel(:,1:ptsnum);
    p2_pixel=p2_pixel(:,1:ptsnum);
    invK=inv(K);
%     p1_noise=[rand(2,ptsnum).*0.1;zeros(1,ptsnum)];
%     p2_noise=[rand(2,ptsnum).*0.1;zeros(1,ptsnum)];
%     p1_noisy=p1_pixel+p1_noise;
%     p2_noisy=p2_pixel+p2_noise;
    [fiveptR,fiveptT]=fivept_RT(p1_pixel,p2_pixel,ptsnum,K);
%     run('find_RT_5ptalgo.m');
    errs(m,:)=compute_reprojection_error_all(fiveptR,fiveptT,p1_pixel,p2_pixel,ptsnum,K);
    fiveptRs{m,1}=fiveptR;
    fiveptTs{m,1}=fiveptT;
end
errs_gt=compute_reprojection_error_all(R_gt,T_gt,p1_pixel,p2_pixel,ptsnum,K);