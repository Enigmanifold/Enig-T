load('GD_epimesh2.mat');
load('sorted_xstd2.mat');
cam1_implane=zeros(4,length(total_cell)-1);
cam2_implane=zeros(4,length(total_cell)-1);
ps_cam1=zeros(2,length(sorted_xstd));
ps_cam2=zeros(2,length(sorted_xstd));
for m=1:length(total_cell)-1
    [x1,y1]=imsphere2implane(total_cell{m+1,1}(1),total_cell{m+1,1}(2));
    [x2,y2]=imsphere2implane(total_cell{m+1,3}(1),total_cell{m+1,3}(2));
    [x3,y3]=imsphere2implane(total_cell{m+1,2}(1),total_cell{m+1,2}(2));
    [x4,y4]=imsphere2implane(total_cell{m+1,4}(1),total_cell{m+1,4}(2));
    cam1_implane(1:2,m)=[x1,y1]';
    cam1_implane(3:4,m)=[x2,y2]';
    cam2_implane(1:2,m)=[x3,y3]';
    cam2_implane(3:4,m)=[x4,y4]';
end
for m=1:length(sorted_xstd)
%     [x_1,y_1]=imsph2implane(sorted_xstd(m,1).*pi./180,sorted_xstd(m,2).*pi./180);
    [x_2,y_2]=imsph2implane(sorted_xstd(m,1).*pi./180,sorted_xstd(m,2).*pi./180);
    ps_cam1(:,m)=[x_1,y_1]';
    ps_cam2(:,m)=[x_2,y_2]';
end
cam1_ori=cam1_implane(1:2,:);
cam1_GD=cam1_implane(3:4,:);
cam2_ori=cam2_implane(1:2,:);
cam2_GD=cam2_implane(3:4,:);
scatter(cam2_ori(1,:),cam2_ori(2,:),'b')
hold on
scatter(cam2_GD(1,:),cam2_GD(2,:),'r')
scatter(e2(1),e2(2),'*','k')
title('mesh points and gradient descent points for Cam2')
xlabel('x')
ylabel('y')
legend('mesh points','gradient descent points','ground truth');
hold off
figure
scatter(ps_cam2(1,:),ps_cam2(2,:),'r')
hold on
scatter(e2(1),e2(2),'*','k')
title('particle swarm results for Cam2')
xlabel('x')
ylabel('y')
legend('particle swarm points','ground truth');
hold off