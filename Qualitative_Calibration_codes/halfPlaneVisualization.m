clear;
clc;
%% Generate random Camera Pose
run('/home/hongyi/Documents/Package/vlfeat/vlfeat/toolbox/vl_setup.m');
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

e1 = null(EGT);
e2 = null(EGT');
e1 = e1 / e1(3);
e2 = e2 / e2(3);

%% Draw Different camera plane
T = T .* 3
c = -rMat' * T;
%Draw Camera center and baseline
plot3(0,0,0,'r.','MarkerSize',20); hold on
plot3(c(1),c(2),c(3),'b.','MarkerSize',20); hold on
counter = 1;
for i = -6:0.1:6
    baseline(:,counter) = (c ./ norm(c)) .* i;
    counter = counter + 1;
end
plot3(baseline(1,:),baseline(2,:),baseline(3,:),'g-','Linewidth',1.5); hold on

%Draw Camera Plane
%Get 4 corner of a camera
imLength = 1;
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
fill3([cc1(1) cc2(1) cc3(1) cc4(1)],...
      [cc1(2) cc2(2) cc3(2) cc4(2)],...
      [cc1(3) cc2(3) cc3(3) cc4(3)],[0.3 0.3 0.3 0.3],'FaceAlpha',.3); hold on

cp1 = rMat' * (cc1 - T);
cp2 = rMat' * (cc2 - T);
cp3 = rMat' * (cc3 - T);
cp4 = rMat' * (cc4 - T);

plot3([cp1(1) cp2(1)],[cp1(2) cp2(2)],[cp1(3) cp2(3)],'b-'); hold on
plot3([cp2(1) cp3(1)],[cp2(2) cp3(2)],[cp2(3) cp3(3)],'b-'); hold on
plot3([cp3(1) cp4(1)],[cp3(2) cp4(2)],[cp3(3) cp4(3)],'b-'); hold on
plot3([cp4(1) cp1(1)],[cp4(2) cp1(2)],[cp4(3) cp1(3)],'b-'); hold on
plot3([c(1) cp1(1)],[c(2),cp1(2)],[c(3),cp1(3)],'b-'); hold on
plot3([c(1) cp2(1)],[c(2),cp2(2)],[c(3),cp2(3)],'b-'); hold on
plot3([c(1) cp3(1)],[c(2),cp3(2)],[c(3),cp3(3)],'b-'); hold on
plot3([c(1) cp4(1)],[c(2),cp4(2)],[c(3),cp4(3)],'b-'); hold on
fill3([cp1(1) cp2(1) cp3(1) cp4(1)],...
      [cp1(2) cp2(2) cp3(2) cp4(2)],...
      [cp1(3) cp2(3) cp3(3) cp4(3)],[0.5 0.5 0.5 0.5],'FaceAlpha',.3); hold on
grid on
axis equal

%Draw Epipoles
plot3(e1(1),e1(2),e1(3),'r*','MarkerSize',20); hold on

e2C1 = rMat' * (e2 - T);
plot3(e2C1(1),e2C1(2),e2C1(3),'b*','MarkerSize',20); hold on
axis equal;

%% Draw Half plane
vec = c; 
vecPrep = [(-vec(2) * 2 - vec(3) * 1)/vec(1); 2; 1];
vecPrep = vecPrep ./ norm(vecPrep);
cA = (c ./ norm(c)) .* -5;
cB = (c ./ norm(c)) .* 5;

vecX = [0 -vec(3) vec(2);
          vec(3)  0 -vec(1);
          -vec(2) vec(1) 0];
RRPrep = expm(vecX * pi / 2);
normal = RRPrep * vecPrep;

epiline1 = cross(normal, [0;0;1]);
epiline2 = cross(normal, rMat' * [0;0;1]);

if dot(epiline1, vecPrep) < 0
    epiline1 = -epiline1;
end
if dot(epiline2, vecPrep) < 0
    epiline2 = -epiline2;
end

h2 = plot3([e1(1), e1(1) + epiline1(1) * 3],...
      [e1(2), e1(2) + epiline1(2) * 3],...
      [e1(3), e1(3) + epiline1(3) * 3], 'm-','Linewidth',5);
 
h3 = plot3([e2C1(1), e2C1(1) + epiline2(1) * 3],...
      [e2C1(2), e2C1(2) + epiline2(2) * 3],...
      [e2C1(3), e2C1(3) + epiline2(3) * 3], 'm-','Linewidth',5);

% plot3([cA(1) cA(1) + normal(1) * 10],...
%       [cA(2) cA(2) + normal(2) * 10],...
%       [cA(3) cA(3) + normal(3) * 10],'m-'); hold on

% plot3([cA(1) cA(1) + vecPrep(1) * 5],...
%       [cA(2) cA(2) + vecPrep(2) * 5],...
%       [cA(3) cA(3) + vecPrep(3) * 5],'c-'); hold on
% 
% plot3([cB(1) cB(1) + vecPrep(1) * 5],...
%       [cB(2) cB(2) + vecPrep(2) * 5],...
%       [cB(3) cB(3) + vecPrep(3) * 5],'c-'); hold on

color = rand();
h = fill3([cA(1) cA(1) + vecPrep(1) * 5 cB(1) + vecPrep(1) * 5 cB(1)],...
      [cA(2) cA(2) + vecPrep(2) * 5 cB(2) + vecPrep(2) * 5 cB(2)],...
      [cA(3) cA(3) + vecPrep(3) * 5 cB(3) + vecPrep(3) * 5 cB(3)],[color color color color],'FaceAlpha',.2); hold on

figure(3)
plot([e1(1) e1(1)+epiline1(1)],...
     [e1(2) e1(2)+epiline1(2)],'m-'); hold on
 
figure(4)
epiline2C2 = rMat * epiline2;
plot([e2(1) e2(1)+epiline2C2(1)],...
     [e2(2) e2(2)+epiline2C2(2)],'m-'); hold on
%
for i = 1:20
    delete(h);
    delete(h2);
    delete(h3);
    vecX = [0 -vec(3) vec(2);
          vec(3)  0 -vec(1);
          -vec(2) vec(1) 0];
%     RR = eye(3) + sin(i * 20) * vecX + (1 - cos(i * 20)) * vecX * vecX;
%     
%     vecPrepA = RR * vecPrep;
    RRR = expm(vecX * i * 0.1);
    
    
    
    vecPrepA = RRR * vecPrep;
    
    vecX = [0 -vec(3) vec(2);
          vec(3)  0 -vec(1);
          -vec(2) vec(1) 0];
RRPrep = expm(vecX * pi / 2);
normal = RRPrep * vecPrepA;
color = rand();
epiline1 = cross(normal, [0;0;1]);
epiline2 = cross(normal, rMat' * [0;0;1]);

if dot(epiline1, vecPrepA) < 0
    epiline1 = -epiline1;
end
if dot(epiline2, vecPrepA) < 0
    epiline2 = -epiline2;
end
color2 = rand(1,3);
figure(1)
h3 = plot3([e1(1), e1(1) + epiline1(1) * 3],...
      [e1(2), e1(2) + epiline1(2) * 3],...
      [e1(3), e1(3) + epiline1(3) * 3], 'color',color2,'Linewidth',5);
 
h2 = plot3([e2C1(1), e2C1(1) + epiline2(1) * 3],...
      [e2C1(2), e2C1(2) + epiline2(2) * 3],...
      [e2C1(3), e2C1(3) + epiline2(3) * 3], 'color',color2,'Linewidth',5);
    h = fill3([cA(1) cA(1) + vecPrepA(1) * 5 cB(1) + vecPrepA(1) * 5 cB(1)],...
          [cA(2) cA(2) + vecPrepA(2) * 5 cB(2) + vecPrepA(2) * 5 cB(2)],...
          [cA(3) cA(3) + vecPrepA(3) * 5 cB(3) + vecPrepA(3) * 5 cB(3)],[color color color color],'FaceAlpha',.2); hold on
    
     figure(3)
     plot([e1(1) e1(1)+epiline1(1)],...
     [e1(2) e1(2)+epiline1(2)],'color',color2); hold on
 
 figure(4)
epiline2C2 = rMat * epiline2;
plot([e2(1) e2(1)+epiline2C2(1)],...
     [e2(2) e2(2)+epiline2C2(2)],'color',color2); hold on
end