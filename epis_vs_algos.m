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
multiepiresults=[];
multiepiresults_total={};
for trial=1:30
    run('generate_RT.m');
    run('setup_episearch.m');
    run('search_epipoles.m');
    run('analyze_epipoles.m');
    episearch_results2(6:9)=episearch_results(5:8);
    episearch_results2(1:4)=episearch_results(1:4);
    episearch_results2(5)=phi;
    multiepiresults(trial,:)=episearch_results2;
    multiepiresults_total{trial,1}=episearch_results2;
    multiepiresults_total{trial,2}=results1;
    multiepiresults_total{trial,3}=results2;
    multiepiresults_total{trial,4}=results3;
    multiepiresults_total{trial,5}=results4;
    disp(num2str(trial));
end
