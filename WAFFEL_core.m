
initp = [gammap;tp;cq0.';cqprime;ETHz;Ebeam];  

if param.progress == 2
    options = odeset('RelTol', param.accuracy,'OutputFcn',@odeplot,'OutputSel',Np+1); %,'Stats','on');   % set solver options
    clf(figure(1))
else
    options = odeset('RelTol', param.accuracy);
end

%% Solve system of equations

zsteps = linspace(0,param.lambdau*param.und_periods,length(param.Kset)+1);
zint = round(range(zsteps(1:2))/(param.lambdau/4));
kk = 0;
for ii = 1:length(zsteps)-1
    tstart = tic;    
    param.K = param.Kset(ii);
    tempsol = ode45(@(z,y) multif_noparax6(z,y,param),[zsteps(ii) zsteps(ii+1)],initp,options);
  
    ztemp = linspace(zsteps(ii),zsteps(ii+1),zint);
    ysol = deval(tempsol,ztemp);
    for jj = 1:length(ztemp)

        kk = kk+1;
        sol.z(kk) = ztemp(jj);
        sol.gamma(kk,1:Np) = ysol(1:Np,jj).';
        sol.t(kk,1:Np) = ysol(Np+1:Np*2,jj).';
        sol.cq(kk,1:param.nfreq) = ysol(2*Np+1:2*Np+param.nfreq,jj);
        sol.THz_energy(kk) = ysol(2*Np+2*param.nfreq+1,jj);    
        sol.bunch_energy(kk) = ysol(2*Np+2*param.nfreq+2,jj);
        
        kKlf = delta_omega/2/pi*(sol.cq(kk,:).'.*exp(1i*k_zq(:)*(ztemp(jj)+param.zoffset))).'*exp(-1i.*omega(:)*(z0+ztemp(jj))/c); 
        sol.THz_peak_power(kk) = c*eps0*A_mode*max(real(kKlf))^2/2;  
        sol.bunch_power(kk) = param.I*sum(sol.gamma(kk,:).*q')/sum(q)*me*c^2/e0;  
        sol.bunching(kk,1:param.nfreq) = sum(exp(-1i*(omega.*sol.t(kk,:).')).*q)/sum(q); 

    end
    
    initp = ysol(:,end);    
    tfinal = toc(tstart);
    if param.progress == 1; fprintf('%.3f sec simulation time for step %i out of %i\n',tfinal,ii,length(zsteps)-1); end
end
