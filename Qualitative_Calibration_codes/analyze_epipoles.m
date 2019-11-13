epi12_error1=-ones(size(results1,1),1);
epi12_error2=-ones(size(results2,1),1);
epi12_error3=-ones(size(results3,1),1);
epi12_error4=-ones(size(results4,1),1);
for m=1:size(results1,1)
    epi1_error=calculate_ang_dist(epi1_sph',results1(m,1:2));
    epi2_error=calculate_ang_dist(epi2_sph',results1(m,3:4));
    epi12_error1(m)=epi1_error+epi2_error;
end
for m=1:size(results2,1)
    epi1_error=calculate_ang_dist(epi1_sph',results2(m,1:2));
    epi2_error=calculate_ang_dist(epi2_sph',results2(m,3:4));
    epi12_error2(m)=epi1_error+epi2_error;
end
for m=1:size(results3,1)
    epi1_error=calculate_ang_dist(epi1_sph',results3(m,1:2));
    epi2_error=calculate_ang_dist(epi2_sph',results4(m,3:4));
    epi12_error3(m)=epi1_error+epi2_error;
end
for m=1:size(results4,1)
    epi1_error=calculate_ang_dist(epi1_sph',results4(m,1:2));
    epi2_error=calculate_ang_dist(epi2_sph',results4(m,3:4));
    epi12_error4(m)=epi1_error+epi2_error;
end
epi12_error1=epi12_error1.*180./pi;
epi12_error2=epi12_error2.*180./pi;
epi12_error3=epi12_error3.*180./pi;
epi12_error4=epi12_error4.*180./pi;
episearch_results=-ones(1,8);
episearch_results(1:4)=epi12_sphgt;
episearch_results(5:8)=[sum(epi12_error1<1)/size(results1,1),sum(epi12_error2<1)/size(results2,1),sum(epi12_error3<1)/size(results3,1),sum(epi12_error4<1)/size(results4,1)];
results1(:,1:4)=results1(:,1:4).*180./pi;
results2(:,1:4)=results2(:,1:4).*180./pi;
results3(:,1:4)=results3(:,1:4).*180./pi;
results4(:,1:4)=results4(:,1:4).*180./pi;