% More advanced algos3D_vs_images.
k1=2;
e1s_sph=setup_2Dmesh(k1);
top_frac=0.06;
PSruns=1;
diff_amplifier=1000000;
options=optimset('MaxFunEvals',10000,'MaxIter',100000,'TolFun',1e-4);
viewSelection=[640,721;150,182;235,247;12,42;327,333;489,510;560,580;194,213;57,67;23,46;512,534;818,842;880,888;466,499;302,322;358,369;415,423;666,688;713,727;813,817];
% viewSelection=[666,688];
impair_num=size(viewSelection,1);
train_ptsnum=60;
test_ptsnum=50;
GDRE_reproj_errs=-ones(impair_num,test_ptsnum);
PSRE_reproj_errs=-ones(impair_num,test_ptsnum);
fivept_reproj_errs=-ones(impair_num,test_ptsnum);
GT_reproj_errs=-ones(impair_num,test_ptsnum);
GDRE_RE_avg=-ones(1,impair_num);
PSRE_RE_avg=-ones(1,impair_num);
fivept_RE_avg=-ones(1,impair_num);
GT_RE_avg=-ones(1,impair_num);
GDRE_reproj_errs2=-ones(impair_num,train_ptsnum);
PSRE_reproj_errs2=-ones(impair_num,train_ptsnum);
fivept_reproj_errs2=-ones(impair_num,train_ptsnum);
GT_reproj_errs2=-ones(impair_num,train_ptsnum);
GDRE_RE_avg2=-ones(1,impair_num);
PSRE_RE_avg2=-ones(1,impair_num);
fivept_RE_avg2=-ones(1,impair_num);
GT_RE_avg2=-ones(1,impair_num);

for m=1:impair_num
    fidCamera = fopen('/Users/tianqitang/CV_Research/images/delivery_area/rig_calibration_undistorted/cameras.txt'); % TODO: Change directory
    fidImages = fopen('/Users/tianqitang/CV_Research/images/delivery_area/rig_calibration_undistorted/images.txt'); % TODO: Change directory
    [viewSel(1),viewSel(2)]=deal(viewSelection(m,1),viewSelection(m,2));
    run('helper.m');
    if K1==K2
        K=K1;
    else
        disp(strcat('K1!=K2 for ',num2str(m)))
        continue
    end
    path1=strcat('/Users/tianqitang/CV_Research/images/delivery_area/images/',images(viewSel(1)).path);
    path2=strcat('/Users/tianqitang/CV_Research/images/delivery_area/images/',images(viewSel(2)).path);
%     imshow(path1)
%     figure
%     imshow(path2);
    image1 = imread(path1);
    image2 = imread(path2);
    image1 = rgb2gray(image1);
    image2 = rgb2gray(image2);
    im1 = single(image1);
    im2 = single(image2);
    if K1==K2
        K=K1;
        invK=inv(K);
    end
    [p1_pixel,p2_pixel,p1,p2,good_corrs]=get_good_corrs(im1,im2,R_gt,T_gt,K1,K2);
    pt1=p1(:,1);
    pt2=p1(:,2);
    pt3=p2(:,1);
    pt4=p2(:,2);
    epi12_sph=setup_3Dmesh2(e1s_sph,ceil(size(e1s_sph,2)/4),pt1,pt2,pt3,pt4);
    valid_pairs=good_corrs(:,1:train_ptsnum);
    test_pairs=good_corrs(:,train_ptsnum+1:train_ptsnum+test_ptsnum);    
    [GDRE_R,GDRE_T,results2]=episearch(epi12_sph,p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,top_frac,PSruns,'GD','RE','pix',diff_amplifier,options,K1,K2);
    [PSRE_R,PSRE_T,results4]=episearch(epi12_sph,p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,top_frac,PSruns,'PS','RE','pix',diff_amplifier,options,K1,K2);
    [fiveptR,fiveptT]=fivept_RT(p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,K1,K2);
    
    GDRE_reproj_errs(m,:)=compute_reprojection_error_all(GDRE_R,GDRE_T,p1_pixel(:,test_pairs),p2_pixel(:,test_pairs),test_ptsnum,K1,K2)*K1(1);
    PSRE_reproj_errs(m,:)=compute_reprojection_error_all(PSRE_R,PSRE_T,p1_pixel(:,test_pairs),p2_pixel(:,test_pairs),test_ptsnum,K1,K2)*K1(1);
    fivept_reproj_errs(m,:)=compute_reprojection_error_all(fiveptR,fiveptT,p1_pixel(:,test_pairs),p2_pixel(:,test_pairs),test_ptsnum,K1,K2)*K1(1);
    GT_reproj_errs(m,:)=compute_reprojection_error_all(R_gt,T_gt,p1_pixel(:,test_pairs),p2_pixel(:,test_pairs),test_ptsnum,K1,K2)*K1(1);
    
    GDRE_reproj_errs2(m,:)=compute_reprojection_error_all(GDRE_R,GDRE_T,p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,K1,K2)*K1(1);
    PSRE_reproj_errs2(m,:)=compute_reprojection_error_all(PSRE_R,PSRE_T,p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,K1,K2)*K1(1);
    fivept_reproj_errs2(m,:)=compute_reprojection_error_all(fiveptR,fiveptT,p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,K1,K2)*K1(1);
    GT_reproj_errs2(m,:)=compute_reprojection_error_all(R_gt,T_gt,p1_pixel(:,valid_pairs),p2_pixel(:,valid_pairs),train_ptsnum,K1,K2)*K1(1);

    GDRE_RE_avg(m)=sum(GDRE_reproj_errs(m,:))/test_ptsnum;
    PSRE_RE_avg(m)=sum(PSRE_reproj_errs(m,:))/test_ptsnum;
    fivept_RE_avg(m)=sum(fivept_reproj_errs(m,:))/test_ptsnum;
    GT_RE_avg(m)=sum(GT_reproj_errs(m,:))/test_ptsnum;
    
    GDRE_RE_avg2(m)=sum(GDRE_reproj_errs2(m,:))/train_ptsnum;
    PSRE_RE_avg2(m)=sum(PSRE_reproj_errs2(m,:))/train_ptsnum;
    fivept_RE_avg2(m)=sum(fivept_reproj_errs2(m,:))/train_ptsnum;
    GT_RE_avg2(m)=sum(GT_reproj_errs2(m,:))/train_ptsnum;
end
figure
plot(linspace(1,impair_num,impair_num),GDRE_RE_avg,'r')
hold on
plot(linspace(1,impair_num,impair_num),PSRE_RE_avg,'g')
plot(linspace(1,impair_num,impair_num),fivept_RE_avg,'b')
plot(linspace(1,impair_num,impair_num),GT_RE_avg,'k')
title(strcat('Algorithm tests on validation sets with',{' '},num2str(train_ptsnum),' training points.'))
xlabel('Trials')
ylabel('Average reprojection errors (pixels)')
legend('Gradient Descent','Particle Swarm','5-point Algorithm','Ground Truth')
hold off

figure
plot(linspace(1,impair_num,impair_num),GDRE_RE_avg2,'r')
hold on
plot(linspace(1,impair_num,impair_num),PSRE_RE_avg2,'g')
plot(linspace(1,impair_num,impair_num),fivept_RE_avg2,'b')
plot(linspace(1,impair_num,impair_num),GT_RE_avg2,'k')
title('Algorithm tests on training sets')
xlabel('Trials')
ylabel('Average reprojection errors (pixels)')
legend('Gradient Descent','Particle Swarm','5-point Algorithm','Ground Truth')
hold off