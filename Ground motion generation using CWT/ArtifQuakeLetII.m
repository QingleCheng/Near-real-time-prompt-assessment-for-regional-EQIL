% ArtifQuakeLetII - CWT based modification of historic records to obtain 
%                   spectrum compatible acceleration series

% Luis A. Montejo (luis.montejo@upr.edu)
% Updates: https://sites.google.com/a/upr.edu/montejo/
% References:
%   Montejo, L. A. and Suarez, L. E. (2013). An improved CWT based algorithm 
% for the generation of spectrum compatible records. 
%   Suarez, L. E., and Montejo, L. A. (2005). Generation of Artificial Earthquakes 
% via the Wavelet Transform. Int. J. Solids and Structures 42, pp. 5905?919.
function a=ArtifQuakeLetII(dsfolder,dsname,eqfolder,eqname,ii)
% clc; clear all; close all; 
a=0;
% required input data:

% dsfolder = 'F:\cql\0-test\ArtifQuakeLetII-Beta\spectra';% directory with the target sectrum
% dsname   = 'Spectra2.txt';                 % name of spectrum file (2 columns
%                                         % periods[s] PSA[g])
% 
% eqfolder = 'F:\cql\0-test\ArtifQuakeLetII-Beta\eqks';   % directory with the accelerograms
% eqname  = 'KMM0061604161026.EW.txt';             % name of earthquake file
outname =['Comp-',ii,'-',eqname];% 'KMM0061604161026.EW-comp.txt'; % output file name

dt   = 0.01;						    % time step of accelerogram [sec]
g    = 1;                               % accel units (if g's input 1);

T1 = 0.1; T2 = 4;         % define period range for matching 
                           % input T1=0 & T2=0 to match the whole spectrum

FF1 = 0.1;
FF2 = 1/(2*dt);        % define the range of frequencies to perform CWT
NS = 100;              % number of scale (s) values to perform CWT

zi  = 0.05;                  % damping ratio for response spectra
nit = 100;                    % max number of iterations

% end of input data
% =========================================================================
% =========================================================================
% load data: 

addpath(dsfolder)           % directory with the design spectra
addpath(eqfolder)           % directory with the seed record

fid_cql=fopen(eqname);
eq = textscan(fid_cql,'%f%f','headerlines', 2);                  % read the design spectrum data file
eq_t=eq{1,1};
eq_a=eq{1,2};
np=length(eq_t);
xg = eq_a'/9.8;                    % copy record in a vector
t  = eq_t';
% eq = load (eqname);                  % read the design spectrum data file
% [nr,nc]  = size(eq); 	             % columns and rows of data file
% np       = nr*nc;				     % original number of data points
% xg(1:np) = eq'/g;                    % copy record in a vector
% np       = length(xg);  
% t  = 0 : dt : (np-1)*dt;

Array=linspace(0,10,1001);
Array_t=Array';
DS = load (dsname);                  % read the design spectrum data file
To = Array_t;%DS(:,1);                        % periods vector (original)
dso = DS(:,1)/9.8;                       % target accels vector (original)
% DS = load (dsname);                  % read the design spectrum data file
% DS = sortrows(DS,1);
% To = DS(:,1);                        % periods vector (original)
% dso = DS(:,2);                       % target accels vector (original)

% wavelet analysis parameters:

omega = pi; zeta  = 0.05;           % wavelet function parameters
freqs  = sort(exp(log(FF1):(log(FF2)-log(FF1))/(NS-1):log(FF2)),'descend');  % frequencies vector
freqs(1) = FF2; freqs(end) = FF1;
T      = 1./freqs;                           % periods vector               
scales = omega./(2*pi*freqs);                % scales vector
nf     = length(freqs);                      % number of frequencies/scales

if (T1<To(1) || T2>To(end)) && (T1~=0&&T2~=0)
    disp('============================================');
    disp('error: period range defined for matching');
    disp('fails outside the target spectrum');
    disp('============================================');
    return
elseif T1==0&&T2==0
    T1 = To(1);
    T2 = To(end);
    if FF2<1/T1
        disp('============================================');
        disp('because of sampling frequency limitations in the seed record')
        disp(['the target spectra can only be matched from T = ',num2str(1/FF2),' s'])
        disp('============================================');
        T1 = 1/FF2;
    end
    if T2>T(end);
        FF1 = 1/T2;
        disp('============================================');
        disp(['FF1 is redefined as ',num2str(FF1),' Hz to match the whole target spectrum']);
        disp('============================================');
        freqs  = sort(exp(log(FF1):(log(FF2)-log(FF1))/(NS-1):log(FF2)),'descend');  % frequencies vector
        freqs(1) = FF2; freqs(end) = FF1;
        T      = 1./freqs;                           % periods vector
        scales = omega./(2*pi*freqs);                % scales vector
        nf     = length(freqs);                      % number of frequencies/scales
    end
end


% Perform Continuous Wavelet Decomposition:

centertime = (ceil(np/2)-1)*dt; % will be used to move the wavelet to the center

centralpart = (ceil(np/2)-1):(ceil(np/2)-1)+np-1; % this vector will be used later to extract
                                                 % the center part of the
                                                 % convolution
                                                 
C = zeros(nf,np); % the CWT coefficients will be stored here

nconv = 2*np-1;                    % size of the conv. signal
nFFTs = pow2(nextpow2(nconv));     % number of pts to calculate the FFTs
Fs    = fft(xg, nFFTs);

for k=1:nf
    wv = zumontw((t-centertime)/scales(k),omega,zeta);
    Fb   = fft(wv, nFFTs);	         
    FsFb = Fs.*Fb;
    convol = ifft(FsFb, 'symmetric')/sqrt(scales(k));
    convol = convol(1:nconv);
    C(k,:) = convol(centralpart);
end

disp('*********')
disp('Wavelet decomposition performed')
disp('*********')

% Generate detail functions:

D    = zeros(nf,np);

for k=1:nf
    wv = zumontw((t-centertime)/scales(k),omega,zeta);
    Fb   = fft(wv, nFFTs);	
    FC   = fft(C(k,:),nFFTs);
    FCFb = FC.*Fb;
    convol = ifft(FCFb, 'symmetric')/(scales(k)^(5/2));
    convol = convol(1:nconv);
    D(k,:) = -convol(centralpart);
end

s = trapz(scales,D); 
ff = max(abs(xg))/max(abs(s));
s  = ff * s;
D  = ff * D;

disp('*********')
disp('Detail functions generated')
disp('*********')

% response spectra from the reconstructed and original signal

PSAxg = responsespectrum(T,xg,zi,dt);
PSAs  = responsespectrum(T,s,zi,dt);

% initial scaling of record:

ds = interp1(To,dso,T)'; % resample target spectrum


if T1>=T(1) && T2<=T(end) && T2>T1
    Tlocs = find(T>=T1&T<=T2);
    sf = sum(ds(Tlocs))./sum(PSAxg(Tlocs)); % initial scaling factor
else
    disp('error in T1 and/or T2 values'); return
end

s  = sf * s; D  = sf * D;

% Iterative Process

meane  = zeros (1,nit+1);
rmse   = zeros (1,nit+1);
hPSAbc = zeros(nf,nit+1);
ns = zeros(np, nit+1);
hPSAbc(:,1) = sf*PSAs;
ns(:,1) = s;
dif = abs( hPSAbc(:,1) - ds ) ./ ds;
meane(1) = mean(dif) * 100;
rmse(1)  = norm(dif) / sqrt(length(dif)) * 100;
factor = ones(1,length(freqs));
DN = D;

for m = 2 : nit+1
    disp(['Now performing iteration ',num2str(m-1),' of ',num2str(nit)])
    factor(Tlocs) = ds(Tlocs)./hPSAbc(Tlocs,m-1) ;
    DN = bsxfun(@times,DN,factor');   
    ns(:,m) = trapz(scales,DN);
    hPSAbc(:,m) = responsespectrum(T,ns(:,m),zi,dt);
    dif = abs( hPSAbc(Tlocs,m) - ds(Tlocs) ) ./ ds(Tlocs);
    meane(m) = mean(dif) * 100;
    rmse(m)  = norm(dif) / sqrt(length(dif)) * 100;
end

[br brloc] = min(rmse);  % locates min error
sc = ns(:,brloc);        % compatible record
PSAsc = hPSAbc(:,brloc); 

% perform baseline correction:

CT = max([1,t(end)/20]); % time to correct
[vel,despl,ccs,cvel,cdespl] = basecorr(t,sc,CT);
kka = 1; 

while sum(isnan(ccs)~=0) ~= 0
    kka = kka + 1;
    CTn = kka*CT;
    if CTn >= t(end/2);
        disp('**baseline correction failed**');
        ccs = sc; cvel=vel; cdespl=despl; break
    end
    [vel,despl,ccs,cvel,cdespl] = basecorr(t,sc,CTn);
end
    
      
PSAccs = responsespectrum(T,ccs,zi,dt);
difin = abs( PSAccs(Tlocs) - ds(Tlocs) ) ./ ds(Tlocs);
meanefin = mean(difin) * 100;
rmsefin  = norm(difin) / sqrt(length(difin)) * 100;

if To(1) == 0; To(1)=0.00999999; end

% figure;subplot(311); plot(t,ccs,'LineWidth',1); grid on;
% axis tight; ylabel('accel. [g]'); title('compatible record')
% subplot(312); plot(t,cvel,'LineWidth',1); grid on;
% axis tight; ylabel('vel./g');
% subplot(313); plot(t,cdespl ,'LineWidth',1); grid on;
% axis tight; xlabel('time [s]'); ylabel('displ./g');
% 
% figure; plot(0:nit,rmse,'-ks','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g',...
%     'MarkerSize',10); grid on; xlabel('Iteration No.');
% ylabel('RMSE [%]'); ylim([0 1.1*rmse(2)]);
% hold on; plot(brloc-1,br,'ro','LineWidth',3,'MarkerSize',12); hold off
% 
% figure ; plot( T,PSAxg, T,PSAccs,'-^', To,dso,'MarkerSize',2 );
% legend(': original record',': compatible record',': target spectrum');grid on
% xlabel('Period [s]'); ylabel('Spectral acceleration [%g]'); 
% xlim([To(1) To(end)]); ylim([0 1.03*max(PSAsc)])
% hold on; plot([T1 T1],[0 1.03*max(PSAsc)],'r--',...
%     [T2 T2],[0 1.03*max(PSAsc)],'r--','linewidth',3); hold off
% 
% figure ; semilogx( 1./T,PSAxg,1./T,PSAccs,'-^', 1./To,dso,'MarkerSize',2 );
% legend(': original record',': compatible record',': target spectrum');grid on
% xlabel('Frequency [Hz]'); ylabel('Spectral acceleration [%g]'); 
% xlim([1/To(end) 1/To(1)]); ylim([0 1.03*max(PSAsc)]) 
% hold on; plot([1/T1 1/T1],[0 1.03*max(PSAsc)],'r--',...
%     [1/T2 1/T2],[0 1.03*max(PSAsc)],'r--','linewidth',3); hold off

QQ = [eqfolder,'\',outname];
save(QQ,'ccs','-ascii')

% figure;subplot(311); plot(t,ccs,'-b',t,zeros(size(t)),'-r','LineWidth',1); 
% axis tight; ylabel('accel. [g]'); title(num2str(outname))
% ylim([-max(abs(ccs)) max(abs(ccs))])
% subplot(312); plot(t,cvel,t,zeros(size(t)),'-r','LineWidth',1); 
% axis tight; ylabel('vel./g');
% ylim([-max(abs(cvel)) max(abs(cvel))])
% subplot(313); plot(t,cdespl,t,zeros(size(t)),'-r','LineWidth',1); 
% axis tight; xlabel('time [s]'); ylabel('displ./g');
% ylim([-max(abs(cdespl)) max(abs(cdespl))])
% 
% figure ; plot( To,dso,'--k','linewidth',3); title(num2str(outname)); hold on;
% plot(T,sf*PSAxg,':r', T,PSAccs,'b-', 'linewidth',2 );
% legend(': target spectrum',': scaled record',': compatible record');
% xlabel('Period [s]'); ylabel('Spectral acceleration [%g]'); 
% xlim([0 T2]); ylim([0 max([PSAsc' sf*PSAxg'])]); hold off
% 
% ao = sf*xg; vo=dt*cumtrapz(ao); do=dt*cumtrapz(vo);
% figure;subplot(311); plot(t,ccs,t,ao,'LineWidth',1); grid on;
% axis tight; ylabel('accel. [g]'); title('original scaled and compatible record')
% subplot(312); plot(t,cvel,t,vo,'LineWidth',1); grid on;
% axis tight; ylabel('vel./g');
% subplot(313); plot(t,cdespl,t,do,'LineWidth',1); grid on; legend(': comp.',': orig.')
% axis tight; xlabel('time [s]'); ylabel('displ./g');

end