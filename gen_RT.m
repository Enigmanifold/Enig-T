function [R,T,e1,e2,phi,K,p1OnImage,p2OnImage,p1ImageOnImage,p2ImageOnImage]=gen_RT
    % clear; clc; close all;
%% Generate random Camera Pose
flag = 0;
while (flag == 0)
% Temp = normrnd(0,1,[3,3]);
% [R,Q] = qr(Temp);
% 
% % R = [0.9418, -0.2457, 0.2296;
% %        0.2602, 0.9649, -0.0344;
% %        -0.2131, 0.0921, 0.9727];
% tP = (rand([1,2]) - 0.5) * 360; %[32,13];
% tP = deg2rad(tP);
% T = [cos(tP(1)) * cos(tP(2)); cos(tP(1)) * sin(tP(2)); sin(tP(1))].*1;%For predefined translation and rotation
% T = [0.630871660447047;0.720028341825237;0.289067699705772]
% T=-T;
%  T = [0.1;0;1];
%  eulRot = [0.2, 0.1, 0.2];
%% Generate random Camera Pose using epipoles and phi.
% u1r=rand*10;
% u2r=rand*10;
% u1th=rand*2*pi;
% u2th=rand*2*pi;
e1_th=rand*pi/2;
e2_th=rand*pi/2;
e1_phi=rand*2*pi;
e2_phi=rand*2*pi;
e1=imsph2implane([e1_th;e1_phi]);
e2=imsph2implane([e2_th;e2_phi]);
% e1=[u1r*cos(u1th);u1r*sin(u1th);1];
% e2=[u2r*cos(u2th);u2r*sin(u2th);1];
phi=rand*2*pi-pi;
T_scale=randsample(2,1)*2-3;
[R1,R2,T1,T2]=epipoles_phi_to_RT(e1,e2,phi,T_scale);
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
imSize = [101*5, 101*5];

% Build Essential Matrix tx * R
tx = [0 -T(3) T(2);
      T(3)  0 -T(1);
      -T(2) T(1) 0];
  
EGT = tx * R;

%%Generate 3D Points 
DX = (rand(1,2000) - 0.5) * 5;
DY = (rand(1,2000) - 0.5) * 5;
DZ = (rand(1,2000) + 1.5) * 2;
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
p1OnImage=p1(:,pointIndex);
p2OnImage=p2(:,pointIndex);
if size(DPointsOnImage,2) > 200
        flag = 1;
    end
end
end