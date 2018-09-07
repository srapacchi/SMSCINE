classdef GT_2DMB_TGRAPPA_leakage < handle & BaseBufferGadget
    
    properties
        
        image_num = 0;
        series_num = 0;
        ksp_size;
        dosave = 0;
        domapT = 0;
		
		NiterSLCgpa = 1;
		kszSLCgpa = [7 7];%[5 4];
		lambdaSLCgpa = 5e-5;%1e-5;
		CalibSz = [32 24];
        time = 1;
    end
    
    methods

        %=========================
        %
        % Config
        %=========================        
        function g = config(g)
            
            fprintf('The resonance frequency is %d\n', g.xml.experimentalConditions.H1resonanceFrequency_Hz);
%             try % this is the only cast from java.lang.Integer that works in Matlab
%                 nc = g.xml.acquisitionSystemInformation.receiverChannels;
%             catch
%                 nc = 1;
%             end
           
            g.image_num = 1;   % todo this needs to be static or global...
            g.series_num = 1;  % todo this needs to be static or global...
            
            % base resolution (final image matrix)
%             nx = g.xml.encoding.encodedSpace.matrixSize.x;
%             ny = g.xml.encoding.encodedSpace.matrixSize.y;
%             nz = g.xml.encoding.encodedSpace.matrixSize.z;

            addpath('~/Documents/MATLAB/ESPIRiT'...
                  ,'~/Documents/MATLAB/ESPIRiT/nufft_files'...
                  ,'~/Documents/MATLAB/ESPIRiT/utils'...
                  ,'~/Documents/MATLAB/ESPIRiT/coilCompression_code'...
                  ,'~/Documents/MATLAB/ESPIRiT/SPIRiT_code');
            
        end
        
        %=========================
        %
        % Process
        %=========================
        function g = process(g, recon_data)

        	Nslices = 1+g.xml.encoding.encodingLimits.slice.maximum;

            imgbuffer =struct(recon_data.data);
	        disp([' Size of buffer :' num2str(size(imgbuffer.data))]);
			refbuffer = struct(recon_data.reference);
	            
	            % store headers
			refheaders = refbuffer.headers;
            imgheaders = imgbuffer.headers;


            
            
			Na= size(imgbuffer.data,5);
            Nb= size(imgbuffer.data,6);

        	if(Nslices>1)

	            disp('Processing MBAND with TGRAPPA (Hadamar)');
	            if(g.time)
	                tic;
	            end
	            if(g.dosave==1)
	                nowstring=datestr(now,'yymmdd_HHMMSS');
	                %targetPath = '~/Documents/MATLAB/stan/data/';
					targetPath = '/tmp/';
	                gxml = g.xml;
	                save(strcat(targetPath,nowstring,'_data'),'gxml','recon_data','-v7.3');
	            end
	            
	            CAIPIshifts = 2*pi*(0:(Nslices-1))/Nslices;

	            % let's start by removing incomplete phases
	            PHSS = [recon_data.data.headers.idx.phase];
	            Nphase = 1;
	            Nlines = sum(PHSS==Nphase);
	            while( sum(PHSS==(Nphase+1))==Nlines )
	            	Nphase = Nphase+1;
	            end
	            %try to pass it on
	            g.xml.encoding.encodingLimits.phase.maximum = double(Nphase);
	            Nphase = Nphase+1; %matlab style


	            % Data
	            recon_data.data.data = recon_data.data.data(:,:,:,:,1:Nphase,:);
	            idxx = [recon_data.data.headers.idx.phase>=Nphase];
	            idxlast = find(idxx,1,'first')-1;
	            if(isempty(idxlast))
	            	idxlast = length(idxx);
	            end
	            mUID = max(recon_data.data.headers.measurement_uid(:));
	            validh = recon_data.data.headers.measurement_uid==mUID;
	            idxfirst = find(validh,1,'first');
	            recon_data.data.headers = recon_data.data.headers.select((idxfirst:idxlast));

	            % Reference
	            mUID = max(recon_data.data.headers.measurement_uid(:));
	            validh = recon_data.reference.headers.measurement_uid==mUID;

	            recon_data.reference.data = recon_data.reference.data(:,:,:,:,1:Nphase,:);
	            idxx = [recon_data.reference.headers.idx.phase>=Nphase];
	            %idxlast = find(idxx,1,'first')-1;
	            recon_data.reference.headers = recon_data.reference.headers.select(logical((1-idxx).*validh));

	            %putBufferQ(g,recon_data.data,recon_data.reference);

				% get GRAPPA acceleration and first sampled ky line
	            GPfact = g.xml.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;

		   		[FirstKy, InitKy, LastKy, FirstKx, LastKx] = FindKspaceBoundaries(imgbuffer.data,g.xml);
		        Nky = length(InitKy:GPfact:LastKy);
	    		Nkx = g.xml.encoding.encodedSpace.matrixSize.x;	
	            
				fprintf(1,'Reduced FOV: InitKy=%d - grappa%d - LastKy=%d\n',InitKy, GPfact,LastKy);
					           
	            % Now use time-grappa for isolating the N slices (using the N first phase for now)

	            %use 7th dimension for slices
	            imgkspace=repmat(recon_data.data.data,[1 1 1 1 1 1 Nslices]);
	            refkspace=repmat(recon_data.reference.data,[1 1 1 1 1 1 Nslices]);
                
                leakkspace=zeros(size(repmat(recon_data.data.data,[1 1 1 1 1 1 Nslices Nslices])));

	            % reload buffers
	            imgbuffer =struct(recon_data.data);
		        disp([' Size of buffer :' num2str(size(imgbuffer.data))]);
				refbuffer = struct(recon_data.reference);
		            
		            % store headers
				refheaders = refbuffer.headers;
	            imgheaders = imgbuffer.headers;

	            
	%          % For some reason, partial Fourier is poorly handled, don't have time to investigate why
	%          kszxml = g.xml.encoding.encodedSpace.matrixSize.y;
	%          dky = kszxml-(LastKy-1)-FirstKy;
	%          imgkspace = circshift(imgkspace,dky,2);


	            for b=1:Nb
	            	% ===============================
	            	% DATA SMS recon (with PAT undersampling)
	            	% ===============================
	            	Nref = floor(Nphase/Nslices)-1;
	            	Ktemp = repmat(recon_data.data.data(:,:,:,:,1),[1 1 1 1 Nref]);
	            	kcalib = repmat(recon_data.data.data(:,:,:,:,1:Nslices,1),[1 1 1 1 1 Nslices]);



		            for s=1:Nslices % slice to separate
		            	for s2=1:Nslices % phases to use for slice to separate
		            		Ktemp = 0*Ktemp; %init
		            		for ref=1:Nref

		            			iphase =(ref-1)*Nslices+s2; 
			            		for s3=iphase:iphase-1+Nslices % phase dimension
			            			fprintf(1,'s=%d /%d- ref=%d / %d- s3=%d /%d\n',s,Nslices, ref,Nref,s3,Nphase);
						            Ktemp(:,:,:,:,ref) = Ktemp(:,:,:,:,ref) + recon_data.data.data(:,:,:,:,s3,b) * exp(1i*2*pi*(s3-1)*(s-1)/Nslices);
						        end
					        end
					        % right now, the reference data is averaged over phases(can be done MUCH better)
					        % kcalib(:,:,:,:,s,s2)=mean(Ktemp,5);

					        %alternative: low rank
					        Ktsz = size(Ktemp);
					        [U,~,~] = svd(reshape(Ktemp,[prod(Ktsz(1:4)) size(Ktemp,5)]),'econ');
					        kcalib(:,:,:,:,s,s2)=reshape(U(:,1),Ktsz(1:4));
				    	end
				    end
				  

				    % kernels = SMS_SPSGcalib(permute(kcalib(FirstKx:LastKx , InitKy:GPfact:LastKy ,:,:,:),[2 1 3 4 5]),g.kszSLCgpa,g.lambdaSLCgpa);
				    startKy = 1+rem(InitKy-1,Nslices);
				    kernels = SMS_SPSGcalib(permute(mycrop(squeeze(kcalib( : , startKy:GPfact:end ,:,:,:,1)),g.CalibSz,'c'),[2 1 3 4 5]),g.kszSLCgpa,g.lambdaSLCgpa);
				    
				    kernels = repmat(kernels,[1 1 1 1 1 Nslices]);
				    for s2 = 2:Nslices
				    	kernels(:,:,:,:,:,s2) = SMS_SPSGcalib(permute(mycrop(squeeze(kcalib( : , startKy:GPfact:end ,:,:,:,s2)),g.CalibSz,'c'),[2 1 3 4 5]),g.kszSLCgpa,g.lambdaSLCgpa);
                    end

                    % RECON using SPLIT SLICE GRAPPA
					for slc=1:Nslices
						for ph=1:Nphase
							for iter=1:g.NiterSLCgpa
		                       imgkspace(FirstKx:LastKx,InitKy:GPfact:LastKy,:,:,ph,b,slc)=...
			                       	 exp(1i*2*pi*(ph-1)*(slc-1)/Nslices)* permute(SMS_SPSGrecon(permute(squeeze(imgkspace(FirstKx:LastKx,InitKy:GPfact:LastKy,:,:,ph,b,slc)),[2 1 3 4])...
			                        , kernels(:,:,:,:,slc,rem(ph-1,Nslices)+1) ),[2 1 3 4]); 
		                    end
		                    %undo CAIPI shift
		                    %imgkspace(:,:,:,:,ph,b,slc) = SMS_CAIPIshift(squeeze(imgkspace(:,:,:,:,ph,b,slc)),CAIPIshifts(slc));

		                end
                    end
                    
                     % EVAL LEAKAGE using SPLIT SLICE GRAPPA
					for slc=1:Nslices
                        for slc2=1:Nslices
                            for ph=1:Nphase
                                %
                                   leakkspace(FirstKx:LastKx,InitKy:GPfact:LastKy,:,:,ph,b,slc,slc2)=...
                                       exp(1i*2*pi*(ph-1)*(slc-1)/Nslices)* permute(SMS_SPSGrecon(permute(squeeze(imgkspace(FirstKx:LastKx,InitKy:GPfact:LastKy,:,:,ph,b,slc)),[2 1 3 4])...
                                        , kernels(:,:,:,:,slc2,rem(ph-1,Nslices)+1) ),[2 1 3 4]); 
 
                                %undo CAIPI shifts
                                leakkspace(:,:,:,:,ph,b,slc,slc2) = SMS_CAIPIshift(squeeze(leakkspace(:,:,:,:,ph,b,slc,slc2)),CAIPIshifts(slc2));

                            end
                        end
                    end
                    

					% ===============================
		            % REF SMS recon (PAT reference)
		            % ===============================
					Ktemp = repmat(recon_data.reference.data(:,:,:,:,1),[1 1 1 1 Nref]);
					kcalib = repmat(recon_data.reference.data(:,:,:,:,1:Nslices,1),[1 1 1 1 1 Nslices]);

		            for s=1:Nslices
	            		Ktemp = 0*Ktemp; %init
	            		for s2=1:Nslices
					        %     Ktemp = Ktemp + recon_data.reference.data(:,:,:,:,s2,b).*exp(1i*2*pi*(s2-1)*(s-1)/Nslices);
					        % end
					        % kcalib(:,:,:,:,s)=Ktemp;
					        for ref=1:Nref
		            			iphase =(ref-1)*Nslices+s2; 
			            		for s3=iphase:iphase-1+Nslices % phase dimension
			            			% fprintf(1,'SSSR REF s=%d /%d- ref=%d /%d- s2=%d /%d\n',s,Nslices, ref,Nref,s2,Nphase);
			            			% disp([' Size of buffer :' num2str(size(Ktemp(:,:,:,:,ref)))]);
			            			% disp([' Size of data   :' num2str(size(recon_data.reference.data(:,:,:,:,s2,b)))]);
						            Ktemp(:,:,:,:,ref) = Ktemp(:,:,:,:,ref) + recon_data.reference.data(:,:,:,:,s3,b) .* exp(1i*2*pi*(s3-1)*(s-1)/Nslices);
						        end
					        end

					        % right now, the reference data is averaged over phases(can be done MUCH better)
					        %kcalib(:,:,:,:,s,s2)=mean(Ktemp,5);
					   		%alternative: low rank
					        Ktsz = size(Ktemp);
					        [U,~,~] = svd(reshape(Ktemp,[prod(Ktsz(1:4)) size(Ktemp,5)]),'econ');
					        kcalib(:,:,:,:,s,s2)=reshape(U(:,1),Ktsz(1:4));
				    	end
				        
				    end

				    if(g.dosave==2)
				    	nowstring=datestr(now,'yymmdd_HHMMSS');
	                	%targetPath = '~/Documents/MATLAB/stan/data/';
	                	targetPath = '/tmp/';
	                	gxml = g.xml;
	                	save(strcat(targetPath,nowstring,'_DEBUG'),'gxml','recon_data','kcalib','Ktemp','-v7.3');
	                end

				    g.CalibSz(2) = min( g.CalibSz(2) , size(kcalib,2) );
				    kernels = SMS_SPSGcalib(permute(mycrop(squeeze(kcalib(:, : ,:,:,:,1)),g.CalibSz,'c'),[2 1 3 4 5]),g.kszSLCgpa,g.lambdaSLCgpa);
				    kernels = repmat(kernels,[1 1 1 1 1 Nslices]);
				    for s2 = 2:Nslices
				    	kernels(:,:,:,:,:,s2) = SMS_SPSGcalib(permute(mycrop(squeeze(kcalib(:, : ,:,:,:,s2)),g.CalibSz,'c'),[2 1 3 4 5]),g.kszSLCgpa,g.lambdaSLCgpa);
				    end
				    %kernels = SMS_SPSGcalib(permute(kcalib(FirstKx:LastKx , InitKy:GPfact:LastKy ,:,:,:),[2 1 3 4 5]),g.kszSLCgpa,g.lambdaSLCgpa);

				    % disp([' Size of ref :' num2str(size(refkspace))]);
				    % disp([' Size of kernels :' num2str(size(kernels))]);
				    % fprintf(1,' kx=%d / %d\n',FirstKx, LastKx);
				    [FirstKy, InitKy, LastKy, FirstKx, LastKx] = FindKspaceBoundaries(refkspace,g.xml);
				    % fprintf(1,' Ref kx=%d / %d\n',FirstKx, LastKx);
					for slc=1:Nslices
						for ph=1:Nphase
							for iter=1:g.NiterSLCgpa
		                       refkspace(FirstKx:LastKx,:,:,:,ph,b,slc)=...
			                       	 exp(1i*2*pi*(ph-1)*(slc-1)/Nslices)* permute(SMS_SPSGrecon(permute(squeeze(refkspace(FirstKx:LastKx,:,:,:,ph,b,slc)),[2 1 3 4])...
			                        , kernels(:,:,:,:,slc,rem(ph-1,Nslices)+1) ),[2 1 3 4]); 
		                    end

		                    %undo CAIPI shift
		                    refkspace(:,:,:,:,ph,b,slc) = SMS_CAIPIshift(squeeze(refkspace(:,:,:,:,ph,b,slc)),CAIPIshifts(slc));
		                end
		            end

		        end % end of b loop


				for slc=1:Nslices
						                

	                %if(slc==1)
	                % modify buffer: change slice number
	                acqtime=imgheaders.acquisition_time_stamp;
            		selected = acqtime>0;
					imgbuffer.headers.idx.slice(logical(selected)) = uint16(slc-1);
            		% remove CINE oversampling ky lines 
            		% oskylines =  imgheaders.idx.kspace_encode_step_1<LastKy;
            		% fprintf('size selected %d oskylines %d \n',size(selected),size(oskylines));
            		% selected = selected.*oskylines;


	               
	                % imgbuffer.headers.position(:,1:Nhdrs)=repmat(refslcpos,[1 Nhdrs]);
	                % TO DO: change slice position!
	                SliceThick = g.xml.encoding.reconSpace.fieldOfView_mm.z; %mm
	                SliceDir = imgbuffer.headers.slice_dir;
	                GUESS_sliceshift  = (slc-1) * 2 * SliceThick * SliceDir; %100 % GAP
	                imgbuffer.headers.position = imgbuffer.headers.position + GUESS_sliceshift;

	                % copy data onto buffer
	                imgbuffer.data = single(imgkspace(:,:,:,:,:,:,slc));

	                Ksz = size(imgbuffer.data);
            		NhdrsNecessary = prod(Ksz(2:end))/Ksz(4);
	                imgbuffer.headers = imgbuffer.headers.select(1:NhdrsNecessary);
	                


	                acqtime=refbuffer.headers.acquisition_time_stamp;
            		selected = acqtime>0;
	                refbuffer.headers.idx.slice(logical(selected)) = uint16(slc-1);
	               
	                % TO DO: change slice position!
	                SliceThick = g.xml.encoding.reconSpace.fieldOfView_mm.z; %mm
	                SliceDir = refbuffer.headers.slice_dir;
	                GUESS_sliceshift  = (slc-1) * 2 * SliceThick * SliceDir; %100 % GAP
	                refbuffer.headers.position = refbuffer.headers.position + GUESS_sliceshift;

	                % copy data onto buffer
	                refbuffer.data = single(refkspace(:,:,:,:,:,:,slc));

	                % let's try by sending the entire reference dataset for now (might need to be cropped to the slice of interest)
	                % pass on data as they are (no processing)
                    
                    g.xml.encodingLimits.set.maximum = uint16(Nslices-1);
                    
                    for slc2=1:Nslices
                        imgbuffer.data = single(leakkspace(:,:,:,:,:,:,slc,slc2));
                        % lets use set dimension for leakage
                        imgbuffer.headers.idx.set(:) = uint16(slc2-1);
                        %refbuffer.headers.idx.set(:) = uint16(slc2-1);
                        
                        putBufferQ(g,imgbuffer,refbuffer);
                        
                    end
                    
                    %end
	            end % end of slice loop

	          
	            if(g.time)
	                fprintf(1,'Slice-GRAPPA time %g\n',toc);
	            end        

	        else %single slice (noMB)
	
	        	% pass on data as they are (no processing)
	        	putBufferQ(g,recon_data.data,recon_data.reference);
	        end


            disp('--- Sending data is done ---');
        
        end %end of process
        
    end %end of methods
end %end of classdef

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
	
	LastKy = min(size(imgkspace,2),LastKy);

end