function [FirstKy, InitKy, LastKy, FirstKx, LastKx] = FindKspaceBoundaries(imgkspace,gxml)

	Kx1=sum(sum(abs(imgkspace(:,:,:,:)),4),1);
    Kx1nnz=find(Kx1);
    FirstKy = Kx1nnz(1);

    Ky1=sum(sum(abs(imgkspace(:,:,:,:)),4),2);
    FirstKx=find(Ky1,1,'first');
    LastKx=find(Ky1,1,'last');
    
    %InitKy = rem(FirstKy,GPfact);
    %if(InitKy==0)
    %    InitKy=GPfact;
    %end
	InitKy=FirstKy;
    LastKy = 1+floor(gxml.encoding.encodedSpace.matrixSize.y/2) ...
        -gxml.encoding.encodingLimits.kspace_encoding_step_1.center...
        +gxml.encoding.encodingLimits.kspace_encoding_step_1.maximum;
    LastKy=find(Kx1,1,'last');
	
	%LastKy = min(size(imgbuffer.data(:,:,1),2),LastKy);

end