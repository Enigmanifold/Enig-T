run('epipole_theta_to_RT.m');
theta2_new=-ones(1,ptsnum);
epsilon=0.001;
check=-ones(1,ptsnum);
t1=T(1)/T(3);
t2=T(2)/T(3);
t3=1;
for m=1:ptsnum
    th2_new=atan(((R(2,2)-t2*R(3,2))*tan(theta1(m))+(R(2,1)-t2*R(3,1)))/((R(1,2)-t1*R(3,2))*tan(theta1(m))+(R(1,1)-t1*R(3,1))));
    if th2_new<0
        th2_new=th2_new+pi;
    end
    theta2_new(m)=th2_new;
    if abs(th2_new-theta2(m))<epsilon || abs(th2_new+pi-theta2(m))<epsilon
        check(m)=1;
    else
        check(m)=0;
    end
end