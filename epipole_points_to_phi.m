run('generate_RT.m');
e1=R'*T;
e1=e1/e1(3);
e2=T/T(3);
ptsnum=size(DPoints,2);
ptsnum=20;
% [phi1,phi2,phi1_alt,phi2_alt]=calc_phi(p1,p2,e1,e2,ptsnum);
phi1=calc_phi2(p1,p2,e1,e2,ptsnum);