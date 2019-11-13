% run('generate_gtepis.m');
% gtepis_size=size(gtepis12_sph,1);
% phi=rand*2*pi;
% algo_results=-ones(gtepis_size,8);
% for m=1:gtepis_size
%     e1=-ones(1,3);
%     e2=-ones(1,3);
%     gtepi12_sph=gtepis12_sph(m,:);
%     [R,T,p1,p2]=episph_phi_to_RT(gtepi12_sph,phi);
%     [e1(1),e1(2)]=imsphere2implane(gtepi12_sph(1),gtepi12_sph(2));
%     [e2(1),e2(2)]=imsphere2implane(gtepi12_sph(3),gtepi12_sph(4));
%     e1(3)=1;
%     e2(3)=1;
%     run('setup_episearch.m');
%     run('search_epipoles.m');
%     run('analyze_epipoles.m');
%     algo_results(m,:)=episearch_results;
% end
trial_num=3;
multiepiresults=[];
multiepiresults_total={};
for trial=1:trial_num
    run('generate_RT.m');
    run('setup_episearch.m');
    run('search_epipoles.m');
    run('analyze_epipoles.m');
    episearch_results2(6:9)=episearch_results(5:8);
    episearch_results2(1:4)=episearch_results(1:4).*180./pi;
    episearch_results2(5)=phi*180/pi;
    multiepiresults(trial,:)=episearch_results2;
    multiepiresults_total{trial,1}=episearch_results2;
    multiepiresults_total{trial,2}=results1;
    multiepiresults_total{trial,3}=results2;
    multiepiresults_total{trial,4}=results3;
    multiepiresults_total{trial,5}=results4;
    disp(num2str(trial));
end
% for trial=1:trial_num
%     block1=multiepiresults_total{trial,1};
%     block1=block1(1:5);
%     block2=multiepiresults_total{trial,2};
%     block3=multiepiresults_total{trial,3};
%     block4=multiepiresults_total{trial,4};
%     block5=multiepiresults_total{trial,5};
%     block1=block1.*180./pi;
%     block2=block2.*180./pi.*180./pi;
%     block3=block3.*180./pi.*180./pi;
%     block4=block4.*180./pi.*180./pi;
%     block5=block5.*180./pi.*180./pi;
%     multiepiresults_total{trial,1}(1:5)=block1;
%     multiepiresults_total{trial,2}=block2;
%     multiepiresults_total{trial,3}=block3;
%     multiepiresults_total{trial,4}=block4;
%     multiepiresults_total{trial,5}=block5;
% end