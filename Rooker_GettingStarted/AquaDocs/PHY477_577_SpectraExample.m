clear all
close all
%% input variables from SIO-pier historical record:
fileIn  = './SIO_pier_historical_temp.mat';
load(fileIn)
% t   =  elapsed time in decimal years
% dt  =  sample interval in decimal years
% Ts  =  detrended surface temperature record (dayly observations)
% N   =  number of observations in record
%
M      = 16;
lap    = 2/3;
[A2 ,f ] = welch_method(Ts,dt,1,0);% function at end-of-file
[A2w,fw] = welch_method(Ts,dt,M,lap);
%
% previously, we estimated that \nu/(2N/Ns)=2
% for overlap=2/3 with Hann-function taper
dof   = 2;
Xp    = chi2inv(0.9,dof);
Xm    = chi2inv(0.1,dof);
err   = (A2w(end)*[dof/Xp  dof/Xm]);
%
dof   = 4*M;
Xp    = chi2inv(0.9,dof);
Xm    = chi2inv(0.1,dof);
err_w = (A2w(end)*[dof/Xp  dof/Xm]);
%
% generate a figure and axes
xm  = 4; ym  = 3; pw  = 20; ph  = 10; ag = 1;
ppos1=[xm ym  pw ph];
cbpos=[xm+0.8*pw ym+0.25*ph ag/2 ph/5];
ps   =[1.5*xm+pw 1.5*ym+ph];
%
% create figure
fig0= figure('units','centimeters','color','w');
pos = get(fig0,'position');
pos(3:4)=ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
f0ax1 = axes('units','centimeters','position',ppos1);
f0ax1p1=loglog(f0ax1,f,A2,'-k',fw,A2w,'-r','linewidth', 2);hold on
f0ax1p2=loglog(f0ax1,f(2)*[1 1], err, 's-k',2*f(2)*[1 1], err_w,'s-r','linewidth',2,'markersize',3);
ylabel(f0ax1,'$S_{\scriptscriptstyle T T}$ [C$^2$/cpy]','interpreter','latex')
xlabel(f0ax1,'$f$ [cycle/year = cpy]','interpreter','latex')
grid on
set(f0ax1,'tickdir','out','ticklabelinterpreter','latex','ylim',[1e-7 1e2],...
          'xlim',[0.5*f(2) f(end)],'xtick',[1e-2 1e-1 1e0 1e1 1e2],'fontsize',25)
%
function [A2,f] = welch_method(x,dt,M,olap);
    % x    = input record of lengh N
    % dt   = sample interval
    % M    = number of sub-windows
    % olap = fractional overalap of sub-windows
    %
    N  = length(x);
    Ns = floor(N/(M*(1-olap)+olap));% Ns elements in window
    ds = floor(Ns-olap*Ns);% shift for start of sub-window
    %
    % reduced frequency resolution
    f0 = 1/(Ns*dt);
    %
    % Positive frequencies are different for even vs. odd N
    if iseven(Ns)
        stop = Ns/2;
    else
        stop = (Ns+1)/2;
    end
    %
    f  = f0*[0:stop-1];
    %
    % each window will be tapered
    win = hanning(Ns);
    % pre-allocate the fft array
    F    = nan(stop,M);
    for m = 1:M
        n   = [1:Ns] + (m-1)*ds;% indices in subwindow
        xp  = x(n);% data in sub-window
        s2  = var(xp);% untapered window variance
        xw  = detrend(xp).*win;% detrend+apply taper
        xw  = xw*sqrt(s2/var(xw));% recover variance
        Xw  = fft(xw);% fast Fourier transform
        F(:,m)  = Xw(1:stop);% log positive frequencies
    end
    % get real/imaginary parts
    C  = real(F);
    Q  = imag(F);
    %
    % average the spectra esimated from each window
    A2 = 2*sum( (C.^2 + Q.^2) ,2)*(dt/(Ns*M));% units=[C^2/cpy]
end
