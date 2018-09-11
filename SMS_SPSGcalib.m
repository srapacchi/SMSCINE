function ws = SMS_SPSGcalib(kdataCal,ksize,lambda)
% Inputs:
%   kdataCal should be ky kx nc mb
%   ksize is the kernels size
%   lambda is the Tykonov regularization (optional)
%
% Output: 
%   ws are the SplitSliceGRAPPA (SPSG) kernels
%
% Performs SPSG calibration
% 

if(~exist('ksize','var') || isempty(ksize) )
    srcx = 3;                         % should be odd
    srcy = 3;    % can now be even and odd
else
    srcx =ksize(1);
    srcy =ksize(2);
end

if(~exist('lambda','var') || isempty(lambda) )
    lambda = 5e-4;
end
[sx,sy,nc,mb] = size(kdataCal);

                     
nyacs=sx; nxacs=sy;

acs_src = kdataCal;
src=zeros(nc*srcy*srcx,(nyacs-(srcy-1))*(nxacs-(srcx-1)));
trg=zeros(nc,(nyacs-(srcy-1))*(nxacs-(srcx-1)));

cnt = 0;
for slcInd = 1:mb
    for xind=floor(srcx./2)+1:nxacs-floor(srcx./2),
        for yind=1:nyacs-(srcy-1),
            cnt=cnt+1;
            src(:,cnt)=reshape(acs_src(yind:yind+(srcy-1),xind-floor(srcx./2):xind+floor(srcx./2),:,slcInd),nc*srcy*srcx,1);
        end
    end
end

% pseudo matrix inversion
% pinv reg
A=src';
AA = A'*A;
S = svd(AA,0);
S = sqrt(max(abs(S)));
X = ((AA+eye(size(AA)).*lambda.*S.^2)\A')';

acs_trg = kdataCal;

ws = zeros([srcy srcx nc nc mb]);

for m=1:mb
    cnt = 0;  % This is a lazy counter. could be done much faster.
    for slcInd = 1:mb
        for xind=floor(srcx./2)+1:nxacs-floor(srcx./2),
            for yind=1:nyacs-(srcy-1),
                cnt=cnt+1;
                if slcInd == m
                    % trg(:,cnt)=reshape(acs_trg(yind+floor(((srcy-1)+1)/2):yind + floor(((srcy-1)+1)/2),xind,:,slcInd),nc,1);
                    trg(:,cnt)=reshape(acs_trg(yind+floor(((srcy-1)+1)/2),xind,:,slcInd),nc,1);
                else
                    trg(:,cnt)=0;
                end
            end
        end
    end
    
    wst = trg*X;
%     ws_(:,:,m) = ws;
    ws(:,:,:,:,m) = reshape(wst.',[srcy srcx nc nc]);
%     for c=1:nc
%         for cc=1:nc
%             tmpres(:,:,cc) = filter2(ws(:,:,cc,c),kdataMB(:,:,cc),'same');
%         end
%         res(:,:,c,m)=sum(tmpres,3);
%     end
    
end



end