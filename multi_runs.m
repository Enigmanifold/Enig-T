runs=100;
errs=-ones(1,runs);
% for xd = 1:runs
    run('epipole_theta_to_RT.m')
    max(error)
%     errs(1,xd)=max(error);
% end