function [avgleak,leakage,mask,MAG] = SMS_calcLeakage(ims)

% function to extract leakage fom GT_2DMB_TGRAPPA_leakage

IMsz= size(ims);

Nslices = IMsz(3);

leakage = zeros(IMsz(1:3));
MAG = zeros(IMsz(1:3));

% make sure time dim (5) is not padded with 0
TimeVals=squeeze(sum(sum(sum(sum(abs(ims),1),2),3),4));
ntrueTime = nnz(TimeVals);
ims=ims(:,:,:,:,1:ntrueTime);

for slc=1:Nslices
    temp = 0*ims(:,:,1,1,:);
    for slc2=1:Nslices
        if(slc2~=slc)
            temp = temp + ims(:,:,slc,slc2,:);
        end
    end 
    MAG(:,:,slc) = median(abs(ims(:,:,slc,slc,:)),5);
    %leakage(:,:,slc)=median(abs(temp),5)./median(abs(ims(:,:,slc,slc,:)),5);
    leakage(:,:,slc)=mean(abs(temp),5)./mean(abs(ims(:,:,slc,slc,:)),5);
    %leakage(:,:,slc)=abs(sum(temp,5))./abs(sum(ims(:,:,slc,slc,:),5));
end

mask = MAG>median(MAG(:));

avgleak = zeros(Nslices,1);
for slc=1:Nslices
    mapI=leakage(:,:,slc);
    
    avgleak(slc)=100*mean(mean(mapI(logical(mask(:,:,slc))),1),2);
end