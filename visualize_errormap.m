clear all;
close all;
run('vary_epi_theta_to_RT.m');
errormap=abs(errormap);
errormapsec=errormap(:,:,1,1);
squeeze_errormapsec=squeeze(errormapsec);
bar3(squeeze_errormapsec);