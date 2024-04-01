%% Use tree kmeans clustering
clear all;close all;
addpath ..\Libraries
addpath ..\DatasetMat
DISP = 0;
REALTIME =1;
S = [];

faceDetector = vision.CascadeObjectDetector();
nose_signal=[];
for i = 1:4
    for j=1:8
        %D = load('BD2-2.mat');
        D = load(sprintf('BD%i-%i.mat',i,j));   
        fprintf(sprintf('BD%i-%i.mat\r\n',i,j));
        D = D.savedata;
        %%
        if(DISP);figure('name','thermalI');end
        k=0;
        D.timeI = D.timeI-D.timeI(1); % change to relative time, index 1 is time =0s
        [M,N]=size(D.ThermalI(:,:,1)); [X,Y] = ndgrid(1:M,1:N);
        maxC = 500; % adaptive color mapping
        traj = [];
        trajVal = [];
        while(k<numel(D.ThermalI(1,1,:)))
            if(k==0); tic(); timenow = 0; end
            
            k=k+1;
            % wait to show it realtime : process frame
            if(DISP && REALTIME)
                while(timenow<D.timeI(k))
                    timenow = toc();
                end
            end
            % read image
           
            %% Create a cascade detector object.
            if(k==1);figure('name','sim','Position',[488 1.8000 560 740.8000]);continue; % skip the first frame
                 I = D.ThermalI(:,:,1); [N,M] = size(I);[X,Y] = ndgrid((1:N)/N,(1:M)/M);  
            else
                    I = D.ThermalI(:,:,k);            In = I-min(I(:)); I = In/max(In(:));
                    %J = integralImage(In);
                    %subplot(3,3,[1 2]); plot(sum(I,1));
                    %subplot(3,3,[6 9]); plot(sum(I,2),numel(sum(I,2)):-1:1);
                    %subplot(3,3,[4 5 7 8]); pcolor(I); shading interp;
                    [id,c]=k_means([5*X(:)';5*Y(:)'; 20*In(:)'],20);
                    pcolor(reshape(id,M,N)); shading interp;
            end

            drawnow;%, hold on, title('Detected features');
             %   plot(points); hold off;
            % Clean up
        end
    end
end