function Kb = SMS_SPSGrecon(Kmb,ws)
% Kmb should be ky kx nc
% ws is the SPSG kernel: nky nkx nc nc
%
% Performs split-slice GRAPPA reconstruction




Kb = 0*Kmb;
Ksz = size(Kmb);
wsz = size(ws);

if(length(Ksz)>3)
    Nim = prod(Ksz(4:end));
    for im=1:Nim
        Kb(:,:,:,im) = SMS_SPSGrecon(Kmb(:,:,:,im),ws);
    end
else
    
    tmpres = 0*Kmb;

    nc=size(Kmb,3);
    for c=1:nc
        for cc=1:nc
            A = filter2(ws(:,:,cc,c),Kmb(:,:,cc),'valid');
            tmpres(ceil(wsz(1)/2):ceil(wsz(1)/2)+size(A,1)-1,ceil(wsz(2)/2):ceil(wsz(2)/2)+size(A,2)-1,cc) = A;
        end
        Kb(:,:,c)=sum(tmpres,3);
    end
end
