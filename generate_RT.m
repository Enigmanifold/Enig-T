clear; clc; close all;
rng(0);
%% Generate random Camera Pose
flag = 0;
while (flag == 0)
Temp = normrnd(0,1,[3,3]);
[R,Q] = qr(Temp);

% R = [0.9418, -0.2457, 0.2296;
%        0.2602, 0.9649, -0.0344;
%        -0.2131, 0.0921, 0.9727];
tP = (rand([1,2]) - 0.5) * 360; %[32,13];
tP = deg2rad(tP);
T = [cos(tP(1)) * cos(tP(2)); cos(tP(1)) * sin(tP(2)); sin(tP(1))].*1;%For predefined translation and rotation
% T = [0.630871660447047;0.720028341825237;0.289067699705772]
% T=-T;
%  T = [0.1;0;1];
%  eulRot = [0.2, 0.1, 0.2];

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