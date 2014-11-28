%% Data analysis using wavelets,
% for github, 28 11 2014
% Make sample M-channel signal, here we have 3 channels
clear all 
Fs=1000;
dt=1/Fs;
f1=8; 
f2=14;
t=(1:dt:2);
s=sin(2*pi.*t*f1);
y=cos(2*pi.*t*f2);

ntrig=3;
nchan=2;
ntimepoints=length(t);


channel1=repmat(s,3,1);
channel2=repmat(y,3,1);
channels=cat(3, channel1, channel2);

data_all=zeros(ntimepoints, nchan, ntrig);
for nc=1:nchan,
    for kk=1:ntrig,
        new1=squeeze(channels(kk,:,nc)); % take channel1 for kk=1
        data_all(:,nc,kk)=new1;
        clear new1;
    end
end

% Evaluate that we have the correct single trials on the correct channel
for nc=1:2,
figure; for kk=1:3, subplot(3,1,kk); plot(squeeze(data_all(:,nc,kk)));end
end



%% Useful params 
ndata=1;
collect_chans_L=data_all;

% The general form  of the input data will be [ntimepoints x nchan x ntrig]
[ntimepoints nchan ntrig]=size(data_all);
% General form of timevector 
%timeVec=(1:length(collect_chans_L)).*1/fs;
timeVec=t;

%% Visualize data
figure;
for cc=1:nchan,
    ctemp=[];
    subplot(nchan, 1, cc); plot(timeVec, data_all(:,cc,1));
    xlabel('time (ms)'); ylabel(['Channel ' num2str(cc)]);
end
suptitle('Our signal - first trial plotted');
clear kk

%% Analysis starts:

%% analysis characteristics %%%%
freqN = input('frequency to start?        ');
repeats = input('For how many times -subsequent frequency bands?   ');
tic
%% Define frequencies
for q = 1:repeats
    freq1=freqN;
    freqN=freq1+20;  
    step=0.2;
    freqVec =freq1:step:freqN; % 2:0.05:16
    disp(freq1)
    disp(freqN)
    width=6;
    for k=1:ndata;
        TFR_array=zeros(length(freqVec), length(timeVec), nchan);
        RHO=0;
        for n=1:nchan  %start for every source - channel
        % Lets select the single trials and detrend them
            collect_sts=zeros(ntimepoints, ntrig-1);   
            collect_sts=data_all(:,n,:);
            collect_sts=squeeze(collect_sts);
            collect_sts=collect_sts'; % ntrig x ntimepoints

            for m=1:ntrig,
                buffer=collect_sts(m, :);
                buffer2=detrend(buffer(:));
                collect_sts(m,:)=buffer2;
                clear buffer buffer2
            end

            B = zeros(length(freqVec), ntimepoints); %% freqVec x timeLength
            PH = zeros(length(freqVec), ntimepoints); %% freqVec x timeLength
            TFR=[]; % empties the variable TFR
            for r=1:ntrig,  % for every single trial     
                for j=1:length(freqVec)  % for every frequency
                    a=squeeze(collect_chans_L(:,n,r))';
                      %[enrg, ph]= energyvec_phase(freqVec(j), collect_sts(r,:),fs, width);
                      %PH(j,:)=ph;
                      %B(j,:)=enrg +B(j,:);
                      B(j, :) = (energyvec(freqVec(j), a, Fs, width)) + B(j,:);
                      clear a
                end % for every frequency
                clear j
                temp(r,n,:,:)=PH(:,:);
            end % for every trigger
            clear r
            TFR = B/ntrig;  % TFR is mean value of B
            
            % or minus one 3marzo2004
            TFR_array(:,:,n) = TFR;
        end   % end for every source -channel
        clear n
  %TFR_all:  ndatasets x nfreqs x ntimepoints x nchan
  TFR_all(k,:,:,:)=TFR_array;
    end % for every data
    clear TFR
    x_array=zeros(ndata, length(timeVec));
    recon_array = zeros(length(freqVec), length(timeVec), nchan);
    for d = 1:nchan
        for f = 1:(length(freqVec))
            x_array = TFR_all(:,f,:,d);  % or TFR_all(k,f,:,d)
            mx_array = mean(x_array, 1);
            recon_array(f,:,d) = mx_array;
        end
        clear f
    end
    clear d
    %end % for repeats
clear q
%% Forming the general name of the dataset 
   %pathname1='C:\Users\Maria\Desktop\1stdatas';
    filename1='sample_data';
    len_1 = length(filename1); 
    save_name=[filename1 '_ANALYSIS_'];    
    stemp1 = [save_name '_' num2str(freq1) '_' num2str(freqN) '_' 'step' num2str(step) '_' 'recon_array_width' num2str(width)]; 
    stempext = ('.mat'); 
    stemp2 = [stemp1 stempext]; 
    pathname_save=pwd; % Or any other path 
    cd(pathname_save)%
    mkdir(stemp1); % 
    cd(stemp1); % 
    eval(['save ' stemp2 ' recon_array freqVec timeVec step filename1 stemp1 s ndata nchan'])
    
    %% Visualize results
    figure;
    for kk=1:nchan,
        subplot(nchan, 1, kk);
        Btemp=recon_array(:,:,kk);
        imagesc(timeVec, freqVec, Btemp);
        xlabel('time (sec)');
        ylabel('frequency (Hz)');
        clear Btemp
    end
    clear kk
    cd(pathname_save)
end  
    %%if we had multiple repeats and datasets here would be the end of
    %%%the first loop
toc    


