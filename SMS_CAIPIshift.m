function K = SMS_CAIPIshift(K,CAIPIshifts,RFphases)
% assumes K is 3D: kx ky coil slice
% CAIPIshifts is a 1D vector of length Nslices (=size(K,4))
% RFphases (optional) is the same dimension as CAIPIshifts (not sure if needed)
%
% Shift image space corresponding to K according to CAIPI shifts
%
% Stan Rapacchi

if(~exist('RFphases','var'))
	RFphases = 0*CAIPIshifts;
end

Ksz = size(K);
Nslices = size(K,4);%Ksz(4);

[xx,~]=meshgrid(1-Ksz(2)/2:Ksz(2)/2,1:Ksz(1));
 xx = repmat(xx,[1 1 Ksz(3)]);
 
for s=1:Nslices
        K(:,:,:,s)=K(:,:,:,s).*exp(1i*CAIPIshifts(s)*xx) * exp(1i*RFphases(s));
end
