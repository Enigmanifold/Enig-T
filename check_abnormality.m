abnormal=[];
for m=1:100
    total_cell_sph=total_total_cell_sph{m};
    if total_cell_sph{1,8}>0.1
        abnormal(end+1)=m;
    end
end