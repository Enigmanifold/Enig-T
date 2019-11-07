% run('setup_episearch.m');
% results=episearch(epi12_sph,p1,p2,ptsnum,top_frac,method,measure,diff_amplifier,options);
epiphi_gt=[epi12_sphgt,best_phi];
epiphi_with_gt=[epiphi;epiphi_gt];
top_frac=0.05;
ptsnum=30;
PSruns=1;
diff_amplifier=1000000;
options=optimset('MaxFunEvals',10000,'MaxIter',100000,'TolFun',1e-4);
% results=calc_phi3_multi(e12s_sph(1,:),p1,p2,ptsnum,diff_amplifier);
% results=episearch(e12s_sph,p1,p2,ptsnum,top_frac,'GD','ReprojectionError',diff_amplifier);
results1=episearch(epi12_sph,p1,p2,ptsnum,top_frac,PSruns,'GD','stdphi',diff_amplifier,options);
results2=episearch(epi12_sph,p1,p2,ptsnum,top_frac,PSruns,'GD','RE',diff_amplifier,options);
results3=episearch(epi12_sph,p1,p2,ptsnum,top_frac,PSruns,'PS','stdphi',diff_amplifier,options);
results4=episearch(epi12_sph,p1,p2,ptsnum,top_frac,PSruns,'PS','RE',diff_amplifier,options);

results_epiphi4=epiphisearch(epiphi_with_gt,p1,p2,ptsnum,top_frac,PSruns,'PS','RE',diff_amplifier,options);