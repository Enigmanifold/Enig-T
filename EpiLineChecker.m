counter = 1;
Greater0 = 0;
Less0 = 0;
for numExp = 1:1000
%% Generate random Camera Pose
flag = 0;
while (flag == 0)
Temp = normrnd(0,1,[3,3]);
[rMat,Q] = qr(Temp);

%rMat = [0.9418, -0.2457, 0.2296;
%        0.2602, 0.9649, -0.0344;
%        -0.2131, 0.0921, 0.9727];
tP = (rand([1,2]) - 0.5) * 360; %[32,13];
tP = deg2rad(tP);
T = [cos(tP(1)) * cos(tP(2)); cos(tP(1)) * sin(tP(2)); sin(tP(1))];%For predefined translation and rotation
%T = [0.7765;0; 2.8978];
%  T = [0.1;0;1];
%  eulRot = [0.2, 0.1, 0.2];

%rMat = eul2rotmatrix(eulRot,'ZYX');
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
  
EGT = tx * rMat;

%%Generate 3D Points 
DX = (rand(1,20000) - 0.5) * 20;
DY = (rand(1,20000) - 0.5) * 20;
DZ = (rand(1,20000) - 0.5) * 20;
DPoints = [DX; 
           DY; 
           DZ];
       
%%Projected on images
for i = 1:size(DPoints,2)
    p1(:,i) = DPoints(:,i) ./ DPoints(3,i);
    p1Image(:,i) = K * p1(:,i);
end

for i = 1:size(DPoints,2)
    P3D2(:,i) = rMat * (DPoints(:,i)) + T;
    p2(:,i) = P3D2(:,i) ./ P3D2(3,i);
    p2Image(:,i) = K * p2(:,i);
end

index1X = p1Image(1,:) < imSize(1) & p1Image(1,:) > 1;
index1Y = p1Image(2,:) < imSize(1) & p1Image(2,:) > 1;
index2X = p2Image(1,:) < imSize(1) & p2Image(1,:) > 1;
index2Y = p2Image(2,:) < imSize(1) & p2Image(2,:) > 1;
index3 = P3D2(3,:) > 0;
index4 = DZ > 0;
pointIndex = index1X & index1Y & index2X & index2Y & index3 & index4;

DPointsOnImage = DPoints(:,pointIndex);
p1ImageOnImage = p1Image(:,pointIndex);
p2ImageOnImage = p2Image(:,pointIndex);

if size(DPointsOnImage,2) > 5
        flag = 1;
    end
end
e11 = EGT(1,1);
e22 = EGT(2,2);
e12 = EGT(1,2);
e21 = EGT(2,1);

e1 = null(EGT);
e2 = null(EGT');

e1 = e1 / e1(3);
e2 = e2 / e2(3);
T = T;


CHECKER = sign(e11 * e22 - e12 * e21);


if CHECKER == -1
    Less0 = Less0 + 1;
else
    Greater0 = Greater0 + 1;
end

figure;
%visualize(rMat, T,e1,e2);
CHECKER

%% Ambiguities 
[R1,T1,R2,T2,R3,T3,R4,T4] = decomposeEssentialMatrix(EGT);



visualize(R1, T1' * 5,e1,e2,'c');
visualize(R2, T2' * 5,e1,e2,'m');
visualize(R3, T3' * 5,e1,e2,'y');
visualize(R4, T4' * 5,e1,e2,'k');


keyboard;


end

function visualize(rMat, T,e1,e2, color)
c = -rMat' * T;
%Draw Camera center and baseline
plot3(0,0,0,'r.','MarkerSize',20); hold on
plot3(c(1),c(2),c(3),'b.','MarkerSize',20); hold on
counter = 1;
for i = -10:0.1:10
    baseline(:,counter) = (c ./ norm(c)) .* i;
    counter = counter + 1;
end
plot3(baseline(1,:),baseline(2,:),baseline(3,:),'g-'); hold on

%Draw Camera Plane
%Get 4 corner of a camera
imLength = 2;
cc1 = [0;0;0] + [0;0;1] + [-imLength;-imLength;0];
cc2 = [0;0;0] + [0;0;1] + [imLength;-imLength;0];
cc3 = [0;0;0] + [0;0;1] + [imLength;imLength;0];
cc4 = [0;0;0] + [0;0;1] + [-imLength;imLength;0];
plot3([cc1(1) cc2(1)],[cc1(2) cc2(2)],[cc1(3) cc2(3)],'r-'); hold on
plot3([cc2(1) cc3(1)],[cc2(2) cc3(2)],[cc2(3) cc3(3)],'r-'); hold on
plot3([cc3(1) cc4(1)],[cc3(2) cc4(2)],[cc3(3) cc4(3)],'r-'); hold on
plot3([cc4(1) cc1(1)],[cc4(2) cc1(2)],[cc4(3) cc1(3)],'r-'); hold on
plot3([0 cc1(1)],[0,cc1(2)],[0,cc1(3)],'r-'); hold on
plot3([0 cc2(1)],[0,cc2(2)],[0,cc1(3)],'r-'); hold on
plot3([0 cc3(1)],[0,cc3(2)],[0,cc1(3)],'r-'); hold on
plot3([0 cc4(1)],[0,cc4(2)],[0,cc1(3)],'r-'); hold on


cp1 = rMat' * (cc1 - T);
cp2 = rMat' * (cc2 - T);
cp3 = rMat' * (cc3 - T);
cp4 = rMat' * (cc4 - T);

plot3([cp1(1) cp2(1)],[cp1(2) cp2(2)],[cp1(3) cp2(3)],color); hold on
plot3([cp2(1) cp3(1)],[cp2(2) cp3(2)],[cp2(3) cp3(3)],color); hold on
plot3([cp3(1) cp4(1)],[cp3(2) cp4(2)],[cp3(3) cp4(3)],color); hold on
plot3([cp4(1) cp1(1)],[cp4(2) cp1(2)],[cp4(3) cp1(3)],color); hold on
plot3([c(1) cp1(1)],[c(2),cp1(2)],[c(3),cp1(3)],color); hold on
plot3([c(1) cp2(1)],[c(2),cp2(2)],[c(3),cp2(3)],color); hold on
plot3([c(1) cp3(1)],[c(2),cp3(2)],[c(3),cp3(3)],color); hold on
plot3([c(1) cp4(1)],[c(2),cp4(2)],[c(3),cp4(3)],color); hold on

%Draw Epipoles
plot3(e1(1),e1(2),e1(3),'r*','MarkerSize',20); hold on

e2C1 = rMat' * (e2 - T);
plot3(e2C1(1),e2C1(2),e2C1(3),'b*','MarkerSize',20); hold on


axis equal
grid on
end