function K = SMS_CAIPIshift(K,CAIPIshifts,RFphases)
% assumes K is 3D: kx ky coil slice
%
% Shift according to CAIPI shifts
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
