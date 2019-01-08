%% now with user input energy distribution
A_e = 2*pi*param.sigmax^2;                                                  % beam cross section 
param.gamma0 = e0*param.Ee/me/c^2;                                                % relativistic gamma factor
deltagamma = param.gamma0*param.deltagammarel;                                    % energy spread 
param.betaz0 = sqrt(1-1/param.gamma0^2-param.Kset(1)^2/2/param.gamma0^2);

if isfield(param,'I')
    param.charge = param.I*param.bunchlength*sqrt(2*pi)/param.betaz0;
else
    param.I = param.charge*param.betaz0/param.bunchlength/sqrt(2*pi);                       % beam current
end
param.ku = 2.*pi./param.lambdau;                                            % undulator wavenumber                                                
param.A_e = A_e;   
param.K = param.Kset(1);


%% initialize phase space (Quiet - start problem )

for ii = 1:2
%     X0 = hammersley(2,param.Np);
%     tp = X0(2,:);  % not actually using any more because of
%     the way q is now assigned

    tp = linspace(0,1,param.Np);
    Np = size(tp',1);

    if (param.gaussian_beam)
        tp = param.bunchlength*8*(tp'-.5);
        %tp = lengthfactor*tp'-lengthfactor/2;
        Nexp = param.I*exp(-tp.^2/2/param.bunchlength^2)*(tp(2)-tp(1))/e0;
        testcharge = param.I*param.bunchlength*sqrt(2*pi)/param.betaz0;
    else
        lengthfactor = param.bunchlength;
        tp = lengthfactor*tp'-lengthfactor/2;
        Nexp = param.I*(tp(2)-tp(1))/e0*ones(Np,1);
    end
    q = poissrnd(Nexp);

    param.Np = ceil(1.15*param.Np^2/length(q(q>0)));
end

%%
dt = min((tp(2)-tp(1))./sqrt(Nexp),lengthfactor);
tau = rand(Np,1).*dt-dt/2;
tp = tp+tau;
param.q =q;

if (param.prebunching)
tp = tp-2.*bunch*param.lambda0/c/2/pi*sin(param.w0*tp+bunchphase);
end

tp(q==0)=[];
q(q==0)=[];
Np = length(tp);

%%
if param.Edist == 0
    gammap = param.gamma0+deltagamma*randn(Np,1);
elseif param.Edist == 1
    gammap = param.gamma0+deltagamma*randn(Np,1);
    betaz = sqrt(1-1./gammap.^2);
    tp = tp-param.tfocus./(betaz*c);
    tp = tp-mean(tp);
elseif param.Edist == 2
    load measEdist
    X0 = hammersley(2,round(Np*length(Epdf)/sum(Epdf)));

    Ecut = interp1(1:length(Epdf),Epdf,linspace(1,length(Epdf),length(X0(1,:))));
    gammap = X0(2,X0(1,:)<Ecut);
    gammap = (gammap-mean(gammap))+param.gamma0;
    chirpeddist = 0;
    if chirpeddist
        if length(gammap) > length(tp)
            cutset = randperm(length(gammap),length(gammap)-length(tp));
            gammap(cutset) = [];
        else
            cutset = randperm(length(tp),length(tp)-length(gammap));
            tp(cutset) = [];    
            q(cutset) = [];
        end
    else
        if length(gammap) > length(tp)
            gammap = gammap(randperm(length(gammap),length(tp)));
        else
            newset = randperm(length(tp),length(gammap));
            tp = tp(newset);    
            q = q(newset);
        end
    end
    betaz = sqrt(1-1./gammap.^2);
    tp = tp.'-param.tfocus./(betaz*c);
    tp = tp.'-mean(tp);
    Np = length(tp);
    gammap = gammap.';
    param.deltagammarel = std(gammap,q)/param.gamma0;
end    

param.Np = Np;
param.q = q;

Ebeam = sum(gammap.*q)*me*c^2;