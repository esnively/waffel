
fig = figure('Name','WafFEL output');
set(fig,'units','normalized','outerposition',[0 0.2 1 .6],'Color',[1,1,1])

axes('Position',[0 0 1 1],'Visible','off');
booltext = {'off','on'};
spectrum = {'Gaussian','EOS'};
parlist = {['Solver accuracy ' num2str(param.accuracy)];' ';
    
    'Radiation Parameters:';['Frozen field approximation ' booltext{param.frozenfield+1}];['# of spectral points ' num2str(param.nfreq)];
    ['Peak frequency ' num2str(c/param.lambda0*1e-12,3) ' THz'];['Bandwidth ' num2str(param.sw*1e-12) ' THz'];['Peak field ' num2str(param.E0*1e-6) ' MV/m'];' ';
    
    'Beam Parameters:';['Space charge ' booltext{param.spacecharge+1}];['# of macroparticles ' num2str(param.Np)];['Charge ' num2str(param.charge*1e12) ' pC']; 
    ['Energy ' num2str(param.Ee*1e-6) ' MeV'];['Energy spread ' num2str(100*param.deltagammarel) ' %'];['Bunchlength ' num2str(param.bunchlength*1e12) ' ps'];' ';
    
    'Undulator and Waveguide:';['Undulator period ' num2str(param.lambdau*1e2) ' cm'];['K parameter ' num2str(param.K)];['# of periods ' num2str(param.und_periods)];
    ['CPPWG spacing ' num2str(param.b*1e3) ' mm']};
text(.01,.5,parlist)

%%

axes('Position',[.12 .6 .12 .34])
[~,peak] = sort(sol.THz_energy);
plot(omega/2/pi,abs(sol.cq(peak(end),1:param.nfreq)).^2,'r','LineWidth',2)
hold on
plot(omega/2/pi,abs(sol.cq(end,1:param.nfreq)).^2,'b:','LineWidth',2)
hold off
xlim([min(omega)/2/pi max(omega)/2/pi])
legend('peak','final','Location','SouthEast')
title('THz Spectrum','FontWeight','bold')
xlabel('frequency (Hz)')

%% Final time-domain field

if param.E0 > 0
    kKlfalt = delta_omega/2/pi*(cq0(:).*exp(1i*k_zq(:)*(param.lambdau*param.und_periods))).'*exp(-1i*omega(:)*(param.lambdau*param.und_periods+z0)/c);
else
    kKlfalt = delta_omega/2/pi*(sol.cq(kk,:).'.*exp(1i*k_zq(:)*(0))).'*exp(-1i.*omega(:)*(z0+ztemp(jj))/c); 
end

axes('Position',[.285 .6 .19 .34])
plot(z0+param.lambdau*param.und_periods,real(kKlf)*1e-6,'b:','LineWidth',2)
title('THz waveform','FontWeight','bold')
ylabel('MV/m')
[Tmaxalt,zrefalt] = max(abs(real(kKlfalt*1e-6)));
[Tmax,zref] = max(abs(real(kKlf*1e-6)));
if max(Tmax,Tmaxalt)>0; ylim([-1.2*max(Tmax,Tmaxalt) 1.2*max(Tmax,Tmaxalt)]); end
ax1 = gca;
ax1.XDir = 'reverse';

%% initial longitudinal phase space
axes('Position',[.52 .6 .21 .34])
scatter(-sol.t(1,1:Np).*1e12,sol.gamma(1,1:Np).*511,50.*q./max(q),-sol.t(1,1:Np).*1e12,'filled')
colorbar
title('LPS initial, head =>','FontWeight','bold')
xlabel('ps')
ylabel('keV')

%% final longitudinal phase space
axes('Position',[.775 .6 .21 .34])
scatter(-sol.t(end,1:Np).*1e12,sol.gamma(end,1:Np).*511,50.*q/max(q),-sol.t(1,1:Np).*1e12,'filled')
colorbar
title('LPS final','FontWeight','bold')
xlabel('ps')
ylabel('keV')

rms_bunchlength0 = sqrt(sum((sol.t(1,1:Np)-mean(sol.t(1,1:Np)))'.^2.*q')./sum(q))*1e12;
rms_bunchlengthf = sqrt(sum((sol.t(end,1:Np)-mean(sol.t(end,1:Np))).^2.*q')./sum(q))*1e12;

% Bunching factor
axes('Position',[.12 .09 .12 .34])
plot(sol.z,abs(sol.bunching(:,ceil(param.nfreq/2))),'b-')
title('Bunching Factor','FontWeight','bold')
xlabel('z (m)')
xlim([min(sol.z) max(sol.z)])
legend(['f_0 = ' num2str(omega(ceil(param.nfreq/2))/2/pi*1e-12,2) ' THz'],'Location','SouthWest')

%Radiation power and energy
axes('Position',[.275 .09 .12 .34])
plot(sol.z,max(sol.gamma.').*511,'r-')
title('Max Bunch Energy','FontWeight','bold')
xlabel('z (m)')
ylabel('keV')

axes('Position',[.435 .09 .12 .34])
plot(sol.z,(sol.bunch_energy-sol.bunch_energy(1)).*10^6,'b-')
title('Bunch Energy (uJ)','FontWeight','bold')
xlabel('z (m)')
xlim([min(sol.z) max(sol.z)])

axes('Position',[.59 .09 .12 .34])
plot(sol.z,(sol.THz_energy)*10^6,'r')
title('THz Energy (uJ)','FontWeight','bold')
xlabel('z (m)')
xlim([min(sol.z) max(sol.z)])

%% energy distribution 
edges=linspace(min(sol.gamma(end,1:Np)*511),max(sol.gamma(end,1:Np)*511),70);
[~,bin]=histc(sol.gamma(end,1:Np)*511,edges);
countf=accumarray(bin',q(:));

edgesp=linspace(min(sol.gamma(1,1:Np)*511),max(sol.gamma(1,1:Np)*511),70);
[~,binp]=histc(sol.gamma(1,1:Np)*511,edgesp);
countp=accumarray(binp',q(:));

if length(param.Kset) == 1 && ~param.tapering
    axes('Position',[.76 .09 .2 .34])
else
    axes('Position',[.735 .09 .1 .34])
    plot(zsteps,param.Kset(1:length(zsteps)),'r-')
    if isfield('param','Kset2')
        hold on
        plot(zsteps,param.Kset2(1:length(zsteps)),'b:')
        hold off
        legend('K manual','K auto')
    end
    if length(param.btaper)>1
        hold on
        plot(zsteps,param.btaper*1e3,'b:')
        hold off
        legend('K','b','Location','SouthWest')
    end
    title('Tapering','FontWeight','bold')
    xlabel('z (m)')
    xlim([min(sol.z) max(sol.z)])
    axes('Position',[.86 .09 .11 .34])
end
plot(edges*1e-3,countf./(edges(2)-edges(1)))
hold on
plot(edgesp*1e-3,countp./(edgesp(2)-edgesp(1)))
xlabel('MeV')
title('Beam Energy Spread')
%{
FWFMfinalE = fullwidthcalc(edges,countf,.1);
FWFMinitialE = fullwidthcalc(edgesp,countp,.1);
legend(['FWFM ' num2str(FWFMinitialE,3) ' keV'],['FWFM ' num2str(FWFMfinalE,3) ' keV'],'Location','NorthEast')

%%
rms_energyspread_fin = std(sol.gamma(end,1:Np)*511);
rms_energyspread_in = std(sol.gamma(1,1:Np)*511);
formatSpec = 'Energy spread rms %.3f keV, FWFM %.3f keV, and bunchlength rms %.3f ps \n';
%fprintf(formatSpec, rms_energyspread_in, FWFMinitialE, rms_bunchlength0);
%fprintf(formatSpec, rms_energyspread_fin, FWFMfinalE, rms_bunchlengthf);

xi = param.K^2 / (4 + 2*param.K^2);                                     % input argument for coupling coefficient
JJ = besselj(0,xi)-besselj(1,xi);                                       % coupling coefficient
lambdar = param.lambda0;
rho1D = (param.I*param.gamma0/(8*pi*IA*A_e)*(param.K*JJ*lambdar/(1+param.K^2/2))^2)^(1/3);
Lgain = param.lambdau/(4*sqrt(3)*pi*rho1D); % assumes vanishing energy spread
Lsat = param.lambdau/rho1D;
Psat = rho1D*param.Ee*param.I;
%}