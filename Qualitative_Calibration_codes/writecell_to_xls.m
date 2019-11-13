for m=1:500
    total_cell_sph=total_total_cell_sph{m};
    writecell(total_cell_sph,strcat('epi_phi(',num2str(m),').xls'));
end