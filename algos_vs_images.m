clc;clear all;close all;
fidCamera = fopen('/Users/tianqitang/CV_Research/images/delivery_area/rig_calibration_undistorted/cameras.txt'); % TODO: Change directory
fidImages = fopen('/Users/tianqitang/CV_Research/images/delivery_area/rig_calibration_undistorted/images.txt'); % TODO: Change directory
viewSelection=[640,721];
run('helper.m');
params = [];
params.ratioTestRatio = 0.75;
params.RANSACIterations = 1000;
params.RANSACInlierDistanceThreshold = 0.5;
path1=strcat('/Users/tianqitang/CV_Research/images/delivery_area/images/',images(viewSelection(1)).path);
path2=strcat('/Users/tianqitang/CV_Research/images/delivery_area/images/',images(viewSelection(2)).path);
% imshow(path1)
% figure
% imshow(path2);
image1 = imread(path1);
image2 = imread(path2);
image1 = rgb2gray(image1);
image2 = rgb2gray(image2);
im1 = single(image1);
im2 = single(image2);
[good_pairs,f1,f2,d1,d2] = find_good_pairs(im1,im2,params);
pts_size=size(good_pairs,2);
p1_pixel_all = [f1(1:2,good_pairs(1,:));ones(1,pts_size)];
p2_pixel_all = [f2(1:2,good_pairs(2,:));ones(1,pts_size)];
load('p1_pixel_all.mat');
load('p2_pixel_all.mat');
p1_all=[(p1_pixel_all(1,:)-ones(1,pts_size).*K1(1,3))./K1(1,1);(p1_pixel_all(2,:)-ones(1,pts_size).*K1(2,3))./K1(2,2);ones(1,pts_size)];
p2_all=[(p2_pixel_all(1,:)-ones(1,pts_size).*K2(1,3))./K2(1,1);(p2_pixel_all(2,:)-ones(1,pts_size).*K2(2,3))./K2(2,2);ones(1,pts_size)];
% p1_all_alt=zeros(3,pts_size);
% p2_all_alt=zeros(3,pts_size);
% invK=inv(K);
% for m=1:pts_size
%     p1_all_alt(:,m)=invK*p1_pixel_all(:,m);
%     p2_all_alt(:,m)=invK*p2_pixel_all(:,m);
% end
valid_pairs=[5,7,8,10,11,12,13,14,15,16,17,18,19,20,22,25,26,28,30,31,39,40,41,42,43,44,45,46,47,48];
test_pairs=[linspace(52,78,78-52+1),linspace(85,103,103-85+1)];
% p1_pixel_all=[p1_all(1,:).*K1(1,1)+ones(1,pts_size).*K1(1,3);p1_all(2,:).*K1(2,2)+ones(1,pts_size).*K1(2,3);ones(1,pts_size)];
% p2_pixel_all=[p2_all(1,:).*K2(1,1)+ones(1,pts_size).*K2(1,3);p1_all(2,:).*K2(2,2)+ones(1,pts_size).*K2(2,3);ones(1,pts_size)];
p1_pixel=p1_pixel_all(:,valid_pairs);
p2_pixel=p2_pixel_all(:,valid_pairs);
p1=p1_all(:,valid_pairs);
p2=p2_all(:,valid_pairs);
ptsnum=length(valid_pairs);
R_gt=R12;
T_gt=T12./norm(T12);
e1=R_gt'*T_gt;
e1=e1/e1(3);
e2=T_gt/T_gt(3);
if K1==K2
    K=K1;
    invK=inv(K);
end
% gtepi1_sph=-ones(2,1);
% gtepi2_sph=-ones(2,1);
% [gtepi1_sph(1),gtepi1_sph(2)]=implane2imsphere(e1);
% [gtepi2_sph(1),gtepi2_sph(2)]=implane2imsphere(e2);
% gtepi12_sph=[gtepi1_sph',gtepi2_sph'];
% for m=1:120
%     figure
%     subplot(1,2,1),imshow(path1)
%     hold on
%     scatter(p1_pixel_all(1,m),p1_pixel_all(2,m),'r','*')
%     hold off
%     subplot(1,2,2),imshow(path2)
%     hold on
%     scatter(p2_pixel_all(1,m),p2_pixel_all(2,m),'r','*')
%     hold off
% end
% for m=1:ptsnum
%     figure
%     subplot(1,2,1),imshow(path1)
%     hold on
%     scatter(p1_pixel(1,m),p1_pixel(2,m),'r','*')
%     hold off
%     subplot(1,2,2),imshow(path2)
%     hold on
%     scatter(p2_pixel(1,m),p2_pixel(2,m),'r','*')
%     hold off
% end
run('setup_episearch.m');
run('search_epipoles.m');
run('find_RT_5ptalgo.m');
epi12_algo1=results1(1,1:4);
epi12_algo2=results2(1,1:4);
epi12_algo3=results3(1,1:4);
epi12_algo4=results4(1,1:4);
[R_algo1,T_algo1,phi_algo1]=epipole_corrs_to_RT(epi12_algo1,p1,p2,ptsnum);
[R_algo2,T_algo2,phi_algo2]=epipole_corrs_to_RT(epi12_algo2,p1,p2,ptsnum);
[R_algo3,T_algo3,phi_algo3]=epipole_corrs_to_RT(epi12_algo3,p1,p2,ptsnum);
[R_algo4,T_algo4,phi_algo4]=epipole_corrs_to_RT(epi12_algo4,p1,p2,ptsnum);
R_errors=[compute_R_diff(R_gt,R_algo1),compute_R_diff(R_gt,R_algo2),compute_R_diff(R_gt,R_algo3),compute_R_diff(R_gt,R_algo4)].*180./pi;
T_errors=[compute_T_diff(T_gt,T_algo1),compute_T_diff(T_gt,T_algo2),compute_T_diff(T_gt,T_algo3),compute_T_diff(T_gt,T_algo4)].*180./pi;
RE_gt=compute_reprojection_error_fl(R_gt,T_gt,p1_all(:,test_pairs),p2_all(:,test_pairs),ptsnum);
RE_1=compute_reprojection_error_fl(R_algo1,T_algo1,p1_all(:,test_pairs),p2_all(:,test_pairs),ptsnum);
RE_2=compute_reprojection_error_fl(R_algo2,T_algo2,p1_all(:,test_pairs),p2_all(:,test_pairs),ptsnum);
RE_3=compute_reprojection_error_fl(R_algo3,T_algo3,p1_all(:,test_pairs),p2_all(:,test_pairs),ptsnum);
RE_4=compute_reprojection_error_fl(R_algo4,T_algo4,p1_all(:,test_pairs),p2_all(:,test_pairs),ptsnum);
RE_fivept=compute_reprojection_error_fl(fiveptR,fiveptT,p1_all(:,test_pairs),p2_all(:,test_pairs),ptsnum);
RE_gt_pix=compute_reprojection_error_all(R_gt,T_gt,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_1_pix=compute_reprojection_error_all(R_algo1,T_algo1,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_2_pix=compute_reprojection_error_all(R_algo2,T_algo2,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_3_pix=compute_reprojection_error_all(R_algo3,T_algo3,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_4_pix=compute_reprojection_error_all(R_algo4,T_algo4,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_fivept_pix=compute_reprojection_error_all(fiveptR,fiveptT,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_gt_norm=norm(RE_gt);
RE_1_norm=norm(RE_1);
RE_2_norm=norm(RE_2);
RE_3_norm=norm(RE_3);
RE_4_norm=norm(RE_4);
RE_fivept_norm=norm(RE_fivept);
RE_gt_norm_pix=norm(RE_gt_pix);
RE_1_norm_pix=norm(RE_1_pix);
RE_2_norm_pix=norm(RE_2_pix);
RE_3_norm_pix=norm(RE_3_pix);
RE_4_norm_pix=norm(RE_4_pix);
RE_fivept_norm_pix=norm(RE_fivept_pix);
REs_norm=[RE_gt_norm,RE_1_norm,RE_2_norm,RE_3_norm,RE_4_norm,RE_fivept_norm];
REs_norm_pix=[RE_gt_norm_pix,RE_1_norm_pix,RE_2_norm_pix,RE_3_norm_pix,RE_4_norm_pix,RE_fivept_norm_pix];