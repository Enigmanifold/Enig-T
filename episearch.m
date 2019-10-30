function results=episearch(epi12_sph,p1,p2,ptsnum,top_frac,PSruns,method,measure,diff_amplifier,options)
    epi12_size=size(epi12_sph,1);
    if isequal(measure,'stdphi')
        if isequal(method,'GD')
            func=@(x)calc_phi3_multi(x,p1,p2,ptsnum,diff_amplifier);
        else
            if isequal(method,'PS')
                func=@(x)calc_phi3(x,p1,p2,ptsnum,diff_amplifier);
            end
        end
    else
        if isequal(measure,'RE')
            if isequal(method,'GD')
                func=@(x)epipole_corrs_to_RE(x,p1,p2,ptsnum,'m');
            else
                if isequal(method,'PS')
                    func=@(x)epipole_corrs_to_RE(x,p1,p2,ptsnum,'s');
                else
                    error('Method Unknown.');
                end
            end
        else
            error('Measure Unknown.');
        end
    end
    nvars=4;
    lb=[0,-pi,0,-pi];
    ub=[pi/2-0.01,pi,pi/2-0.01,pi];
    epi12_sph_ordered=-ones(epi12_size,5);
    if and(isequal(method,'GD'),isequal(measure,'stdphi'))
        for m=1:epi12_size
            x=epi12_sph(m,:);
            min_val_initial=norm(func(x));
            epi12_sph_ordered(m,:)=[x,min_val_initial];
        end
        epi12_sph_ordered=sortrows(epi12_sph_ordered,5);
        epi12_sph_ordered=epi12_sph_ordered(1:round(epi12_size*top_frac),:);
        epi12_size2=size(epi12_sph_ordered,1);
        results=-ones(epi12_size2,5);
        for m=1:epi12_size2
            x=epi12_sph_ordered(m,1:4);
            [xopt,fopt]=lsqnonlin(func,x,lb,ub,options);
            results(m,:)=[xopt,fopt];
        end
        results=sortrows(results,5);
    end
    if and(isequal(method,'GD'),isequal(measure,'RE'))
        for m=1:epi12_size
            x=epi12_sph(m,:);
            min_val_initial=norm(func(x));
            epi12_sph_ordered(m,:)=[x,min_val_initial];
        end
        epi12_sph_ordered=sortrows(epi12_sph_ordered,5);
        epi12_sph_ordered=epi12_sph_ordered(1:round(epi12_size*top_frac),:);
        epi12_size2=size(epi12_sph_ordered,1);
        results=-ones(epi12_size2,5);
        for m=1:epi12_size2
            x=epi12_sph_ordered(m,1:4);
            [xopt,fopt]=lsqnonlin(func,x,lb,ub);
            results(m,:)=[xopt,fopt];
        end
        results=sortrows(results,5);
    end
    if and(isequal(method,'PS'),isequal(measure,'stdphi'))
        for m=1:epi12_size
            x=epi12_sph(m,:);
            min_val_initial=func(x);
            epi12_sph_ordered(m,:)=[x,min_val_initial];
        end
        epi12_sph_ordered=sortrows(epi12_sph_ordered,5);
        epi12_sph_ordered=epi12_sph_ordered(1:round(epi12_size*top_frac),:);
        epi12_size2=size(epi12_sph_ordered,1);
        options = optimoptions('particleswarm','InitialSwarmMatrix',epi12_sph_ordered(:,1:4),'SwarmSize',epi12_size2);
        results=-ones(PSruns,5);
        for m=1:PSruns
            [results(m,1:4),results(m,5)]=particleswarm(func,nvars,lb,ub,options);   
        end
        results=sortrows(results,5);
    end
    if and(isequal(method,'PS'),isequal(measure,'RE'))
        for m=1:epi12_size
            x=epi12_sph(m,:);
            min_val_initial=func(x);
            epi12_sph_ordered(m,:)=[x,min_val_initial];
        end
        epi12_sph_ordered=sortrows(epi12_sph_ordered,5);
        epi12_sph_ordered=epi12_sph_ordered(1:round(epi12_size*top_frac),:);
        epi12_size2=size(epi12_sph_ordered,1);
        options = optimoptions('particleswarm','InitialSwarmMatrix',epi12_sph_ordered(:,1:4),'SwarmSize',epi12_size2);
        results=-ones(PSruns,5);
        for m=1:PSruns
            [results(m,1:4),results(m,5)]=particleswarm(func,nvars,lb,ub,options);   
        end    
        results=sortrows(results,5);
    end
end