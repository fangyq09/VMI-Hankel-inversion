function out = VMIHT(ims,ctX,ctY,SL)

Xcoo = -SL:SL;
Ycoo = -SL:SL;
Tmax = 0.05;
DT = 0.001*Tmax;
T = -Tmax:DT:Tmax;
Tp = 0:DT:Tmax;
LenTp = length(Tp);
LenT = length(T);
HXcoo = Xcoo;
MAT = zeros(LenTp,length(HXcoo));
MAT = 2*pi.*transpose(Tp)*HXcoo;
MAT = besselj(0,MAT);
MAT =diag(Tp)*MAT;

tau = 5;
%kc = 1/27060; %%1/(0.25*658)^2
kc = 1/27600; %%xenon calibration
Nmax = floor(SL.^2./(tau.^2));
R = tau.*sqrt(0:Nmax);
E = kc.*(R.^2);

DTheta = pi./360;
Thetas = 0:DTheta:pi;
Xp = transpose(R)*cos(Thetas);
Yp = transpose(R)*sin(Thetas);

NoF = size(ims,3);
clmin = ctX - SL;
clmax = ctX + SL;
romin = ctY - SL;
romax = ctY + SL;

for ii = 1:NoF
    IMDAT = ims(romin:romax,clmin:clmax,ii);
    %raw_ims(:,:,ii) = IMDAT;
    % do Fourier transform
    FT = nufft(transpose(IMDAT),Xcoo,T);
    FT = real(FT);
    FT = transpose(FT);
    % do Hankel transform
    HT = FT(:,LenT-LenTp+1:LenT)*MAT;
		HT = DT.*HT;
    % the intensity should be a nonnegative function
    HT(HT<0) = 0;
    hinv(:,:,ii) = HT;
    % energy spectra
    HT = rot90(HT);
    [grid_X,grid_Y] = meshgrid(Ycoo,HXcoo);
    Iv = interp2(grid_X,grid_Y,HT,Xp,Yp);
    % Iv(isnan(Iv)) =0;
    Iv = Iv*diag(sin(Thetas));
    speed = sum(transpose(Iv),'omitnan');
    speed = DTheta.*R.*speed;
		IE(:,ii) = speed;
end
out = struct('E',E,'IE',IE,'hinv',hinv);
end
