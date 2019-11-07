clear all;
expt_num=100;
reproj_errs_norm=-ones(1,expt_num);
for m=1:expt_num
    run('generate_RT');
    e1_sph=-ones(2,1);
    e2_sph=-ones(2,1);
    [e1_sph(1),e1_sph(2)]=implane2imsphere(e1);
    [e2_sph(1),e2_sph(2)]=implane2imsphere(e2);
    e12phi=[e1_sph',e2_sph',phi];
    ptsnum=100;
    reproj_err_norm=max(epiphi_to_RE(e12phi,p1,p2,ptsnum));
    reproj_errs_norm(m)=reproj_err_norm;
end
max(reproj_errs_norm)