function [R,T]=epipoles_phi_to_RT(e1,e2,phi,R2flip_threshold,T_scale)
%% Construct R2.
v1=e1/norm(e1);
v2=e2/norm(e2);
u2=cross(v1,v2);
rot1=v1'*v2;
u2_cross=[0,-u2(3),u2(2);u2(3),0,-u2(1);-u2(2),u2(1),0];
R2=eye(3)+u2_cross+u2_cross*u2_cross/(1+rot1);
rndnum=rand;
if rndnum>R2flip_threshold
    R2_additional_axis_y=-(e2(1)+1)/e2(2);
    R2_additional_cross=[0,-1,R2_additional_axis_y;1,0,-1;-R2_additional_axis_y,1,0];
    R2_additional_cross2=R2_additional_cross*R2_additional_cross;
    R2_additional=eye(3)+2*e2(2)^2/(2*e2(2)^2+e2(1)^2+2*e2(1)+1)*R2_additional_cross2;
    R2=R2_additional*R2;
end
%% Construct R3.
v2_cross=[0,-v2(3),v2(2);v2(3),0,-v2(1);-v2(2),v2(1),0];
v2_cross2=v2_cross*v2_cross;
R3=eye(3)+sin(phi)*v2_cross+(1-cos(phi))*v2_cross2;
%% Construct R.
R=R3*R2;
%% Construct T.
T=T_scale.*v2;