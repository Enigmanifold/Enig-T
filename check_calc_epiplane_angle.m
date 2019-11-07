run('generate_RT');
total_pts=size(p1,2);
epiplane_angles1=-ones(1,total_pts-1);
epiplane_angles2=-ones(1,total_pts-1);
e1_sph=zeros(2,1);
e2_sph=zeros(2,1);
[e1_sph(1),e1_sph(2)]=implane2imsphere(e1);
[e2_sph(1),e2_sph(2)]=implane2imsphere(e2);
for m=1:total_pts-1
    pt1_cam1=p1(:,1);
    pt1_cam2=p2(:,1);
    pt2_cam1=p1(:,m+1);
    pt2_cam2=p2(:,m+1);
    epiplane_angles1(m)=calc_epiplane_angle(e1_sph,pt1_cam1,pt2_cam1);
    epiplane_angles2(m)=calc_epiplane_angle(e2_sph,pt1_cam2,pt2_cam2);
end
inds=find(abs(epiplane_angles1-epiplane_angles2)>1e-10);
error_angle_sum=epiplane_angles1(inds)+epiplane_angles2(inds);