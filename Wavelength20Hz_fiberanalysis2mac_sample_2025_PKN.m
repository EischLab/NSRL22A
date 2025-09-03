   
%%Import Fiber Photometry File, delete the first row, and Save as Array A. Also copy the filename to basename as a string. 
% Fiber photometry data and video data saved in .txt and .csv., respectively. 
% For details, ask Sanghee.  
%
% open fiber data and store as number array A.
% open video data and store as cell array B. (Cell array handles mixed data
% types like number and string). 
%
% A: Fiber data (cleaned up raw)
% AA: Fiber data (temp raw data) 
% B: video data (raw)
% C: Fiber + video data (processed, full)
% CC: Fiber + video data (processed, short) 
% D: Fiber (z-score converted, just for peak detection)
% E: video data (cleaned,processed) 
% 
% Any question? contact Haji - takanoh@chop.edu  
%
close all
clear all

%select a fiberphotmetry file and read.
[filename,pathname]=uigetfile('*.txt');
addpath pathname;
cd(pathname);
openfilename=fullfile(pathname,filename);
fiberDataRaw=dlmread_empty(openfilename,'\t',1,0,NaN);

basename1=strrep(filename,'.txt','');

lengthFiberData=size(fiberDataRaw,1)

% clean up array for NaN (empty array)
counter=1;
for iii=1:lengthFiberData;
   nanornot=sum(fiberDataRaw(iii,:));
   if isnan(nanornot)== 0;
       fiberData(counter,:)=fiberDataRaw(iii,:);
       counter=counter+1;
   else
        fiberData=fiberDataRaw;
   end
end

lengthFiberData=size(fiberData,1);
%%
%select a video data
[filename2,pathname2]=uigetfile('*.csv');
basename2=strrep(filename2,'.csv','');
addpath pathname2;1
cd(pathname2);
openfilename=fullfile(pathname2,filename2);


% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = [" ", ","];

% Specify column names and types
opts.VariableNames = ["Datapoint", "Frames", "X", "Y", "Corner1", "Corner2", "InteractionZone"];
opts.VariableTypes = ["double", "double", "double", "double", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, ["Corner1", "Corner2", "InteractionZone"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Corner1", "Corner2", "InteractionZone"], "EmptyFieldRule", "auto");

% Import the data
Btemp = readtable(openfilename, opts);

% Convert to output type
videoDataRaw = table2cell(Btemp);
numIdx = cellfun(@(x) ~isnan(str2double(x)), videoDataRaw);
videoDataRaw(numIdx) = cellfun(@(x) {str2double(x)}, videoDataRaw(numIdx));

clear opts
%%
NSegment=1;
startpoint=1;
endpoint=lengthFiberData;
    
    vlength=size(videoDataRaw,1);
    videoDataClean=zeros(vlength,15);%% E1=1,2,3, E2=0.033,0.066 in sec, E3=original zone binary , E4 reversed binary E5,processed 
    videoDataClean(:,1)=[1:vlength];
    videoDataClean(:,2)=0.033333*(videoDataClean(:,1)-1);
    videoDataClean(:,4)=1;
    videoDataClean(:,8)=1;
    videoDataClean(:,12)=1;
    
% 1st Zone --> 1st Zone     
% Identify       
    for iii=1:vlength;
        if videoDataRaw{iii,5}=="True";
             videoDataClean(iii,3)=1; % zone
             videoDataClean(iii,4)=0; % reverse (non zone) for cleaning 
        elseif videoDataRaw{iii,5}=="TRUE";
             videoDataClean(iii,3)=1; % zone
             videoDataClean(iii,4)=0; % reverse (non zone) for cleaning 
        else
        end
    end
       
% cleaning up interaction zone  - short duration break from interacting will be ignored and filled as interacting     
    videoDataClean(:,4)=CountConsecutiveOnes(videoDataClean(:,4)); % custom consecutiveone function 
    for iii=1:vlength;
        if videoDataClean(iii,4)<10; % 330ms 
            videoDataClean(iii,4)=0;
        else
        end
    end
    
% E(:,5) is reverse of reverse     
    for iii=1:vlength;
        if videoDataClean(iii,4) == 0;
            videoDataClean(iii,5)=1;
        else
        end
    end
 % then, cleanup one more time to remove short duration interaction    
    videoDataClean(:,5)=CountConsecutiveOnes(videoDataClean(:,5));   
    for iii=1:vlength;
        if videoDataClean(iii,5)<15;   % 500ms 
            videoDataClean(iii,5)=0;
        else
        end
    end
    
  % simply binarize the cleanup results one more time in a new column   
    for iii=1:vlength;
        if videoDataClean(iii,5)>0;
            videoDataClean(iii,6)=1;
        else
        end
    end
    
  
 % 2nd Zone --> 2nd Zone    
  % identify Corner 1 zone     
    for iii=1:vlength;
        if videoDataRaw{iii,6}=="True";
             videoDataClean(iii,7)=1; % zone
             videoDataClean(iii,8)=0; % reverse (non zone) for cleaning 
        elseif videoDataRaw{iii,6}=="TRUE";
             videoDataClean(iii,7)=1; % zone
             videoDataClean(iii,8)=0; % reverse (non zone) for cleaning 
        else
        end
    end
    
    % cleaning up corner zone1  - short duration break from interacting will be ignored and filled as interacting     
    videoDataClean(:,8)=CountConsecutiveOnes(videoDataClean(:,8)); % custom consecutiveone function 
    for iii=1:vlength;
        if videoDataClean(iii,8)<10; % 330ms 
            videoDataClean(iii,8)=0;
        else
        end
    end
    
% E(:,5) is reverse of reverse     
    for iii=1:vlength;
        if videoDataClean(iii,8) == 0;
            videoDataClean(iii,9)=1;
        else
        end
    end
 % then, cleanup one more time to remove short duration interaction    
    videoDataClean(:,9)=CountConsecutiveOnes(videoDataClean(:,9));   
    for iii=1:vlength;
        if videoDataClean(iii,9)<15;   % 500ms 
            videoDataClean(iii,9)=0;
        else
        end
    end
    
  % simply binarize the cleanup results one more time in a new column   
    for iii=1:vlength;
        if videoDataClean(iii,9)>0;
            videoDataClean(iii,10)=1;
        else
        end
    end
    
    
    
    
 % 3rd Zone --> 3rd Zone   
    % identify interaction zone     
    for iii=1:vlength;
        if videoDataRaw{iii,7}=="True";
             videoDataClean(iii,11)=1; % zone
             videoDataClean(iii,12)=0; % reverse (non zone) for cleaning 
        elseif videoDataRaw{iii,7}=="TRUE";
             videoDataClean(iii,11)=1; % zone
             videoDataClean(iii,12)=0; % reverse (non zone) for cleaning 
        else
        end
    end
    
    % cleaning up interaction zone  - short duration break from interacting will be ignored and filled as interacting     
    videoDataClean(:,12)=CountConsecutiveOnes(videoDataClean(:,12)); % custom consecutiveone function 
    for iii=1:vlength;
        if videoDataClean(iii,12)<10; % 330ms 
            videoDataClean(iii,12)=0;
        else
        end
    end
    
% E(:,5) is reverse of reverse     
    for iii=1:vlength;
        if videoDataClean(iii,12) == 0;
            videoDataClean(iii,13)=1;
        else
        end
    end
 % then, cleanup one more time to remove short duration interaction    
    videoDataClean(:,13)=CountConsecutiveOnes(videoDataClean(:,13));   
    for iii=1:vlength;
        if videoDataClean(iii,13)<15;   % 500ms 
            videoDataClean(iii,13)=0;
        else
        end
    end
    
  % simply binarize the cleanup results one more time in a new column   
    for iii=1:vlength;
        if videoDataClean(iii,13)>0;
            videoDataClean(iii,14)=1;
        else
        end
    end
    
    
    
    
    
    
%    %corner zone 
%    for iii=1:vlength;
%         if videoDataRaw{iii,5}=="True" | videoDataRaw{iii,6}=="True";
%              videoDataClean(iii,7)=1; % corner zones
%          elseif videoDataRaw{iii,5}=="TRUE" | videoDataRaw{iii,6}=="TRUE";
%              videoDataClean(iii,7)=1; % corner zones
%         end
%    end
   %other area
   for iii=1:vlength;
        if videoDataClean(iii,6)==0 && videoDataClean(iii,10)==0 && videoDataClean(iii,14)==0 ;
             videoDataClean(iii,15)=1; % other area
         else
        end
    end
    
 %   

 %% Run for different control for peak detection +  
 % Run to create full figure analysis

 reply1 = input(['Do you want normal peak detection using original 415nm data (1),' ...
     'using smoothed 415 nm data (2), or using only 470nm data (3):'],"s");
  
    
 % Fiberphotometry data (fluorescecne@470nm)  - baseline correction by fluorescecne@410nm.  
 % initilize array
 
    Fitted=[];
    DeltaFoverF=[];
    combinedFibVidData=zeros(lengthFiberData,11);% output array % edited 11/10/23, 11th colummn 
    fiberZScores=zeros(lengthFiberData,2);
    finalCombinedData=zeros(lengthFiberData,8); % edited 11/10/23 added 8th columnn 
    
    %Control array is 415nm data. Raw array is 470nm data. 
    tempControl=fiberData(:,4);
    Raw=fiberData(:,5);

   if isempty(reply1)
      reply1 = '1';
   end
   
   if reply1 =='1'
       Control = smoothdata(tempControl, 'rlowess', 5);
   elseif reply1 == '2'
       Control = smoothdata(tempControl, 'rlowess', 500);
   elseif reply1 == '3'
       Control = sgolayfilt(Raw, 3, 501);
%        Control = smoothdata(Raw, 'rlowess', 300);
   end
    
    H0=figure('Name',basename1);
    figure(H0);
    subplot(5,2,1)
    plot(tempControl,'Color','blue');
    hold on
    plot(Raw,'Color','green');
    hold on
    plot(Control,'Color','magenta');
    
    combinedFibVidData(:,1:5)=fiberData(:,1:5);
    fiberZScores(:,1)=fiberData(:,1);
    combinedFibVidData(:,3)=0.05*(combinedFibVidData(:,1)-1);
    finalCombinedData(:,1:3)=combinedFibVidData(:,1:3);
    
    %Linear fit model 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     IniValues=[1.055 0.0012];
%     options = optimset('Display','notify','MaxIter',100000,'MaxFunEvals',100000);
%     Estimates = fminsearch(@myfitFPSY,IniValues,options,Control, Raw);
%     %%%Results 
%     a=Estimates(1);
%     b=Estimates(2);
% 
%     %%%Subtract fitted 410 from 470
%     Fitted(:,1)=(a*Control+b);
%     FittedMean=mean(Fitted(:,1));% use this single number when dividing 
%     DeltaFoverF=(Raw-Fitted(:,1))/FittedMean; % This is DF/F0 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
  
    
    SegmentLength=floor(lengthFiberData/NSegment);
    SubtractedA=zeros(lengthFiberData,1);%% Preallocate
    FittedSegment=zeros(lengthFiberData,1);%% Preallocate
%     SegmentFrames=zeros(NSegment,1);
    
    SegmentPoints=ones(NSegment,1);

    for ii=1:NSegment;
    data415=[];
    data470=[];
    
        if ii ~= NSegment;
    mm = 1+SegmentLength*(ii-1);
    nn = SegmentLength*(ii);
    SegmentPoints(ii,1)=mm;
         else
    mm = 1+SegmentLength*(ii-1);
    nn= lengthFiberData;
    SegmentPoints(ii,1)=mm;
        end
        
    data415=Control(mm:nn);
    data470=Raw(mm:nn);

%     IniValues=[1.00 -0.002];
%     options = optimset('Display','notify','MaxIter',100000,'MaxFunEvals',100000);
% %     Estimates = fminsearch(@myfitFPSY,IniValues,options,TempC, TempA);
%     Estimates = fminsearchbnd(@myfitFPSY,IniValues,[0.100 -inf],[5.000 inf],options,TempC, TempA);
%     %%%Results 
%     a=Estimates(1);
%     b=Estimates(2);

    
   
    
    P = polyfit(data415, data470, 3);

    %%%Subtract Iso from Calcium 
    TempFitted=P(1)*data415.^3 + P(2)*data415.^2 + P(3)*data415 + P(4);
    TempSubtracted=(data470 - TempFitted);
    
    SubtractedA(mm:nn)=TempSubtracted;
    FittedSegment(mm:nn)=TempFitted;
  


            
    end
    figure(H0);
    subplot(5,2,[3,4])
    plot(FittedSegment,'LineWidth',2, 'Color','green');
    hold on;
    plot(Raw,'Color',[0.4 0.4 0.4]);
    
    hold on
    for kk=1:NSegment;
     XX=[SegmentPoints(kk,1)  SegmentPoints(kk,1)];
     YY=[min(Raw)  max(Raw)];
     plot(XX,YY,'Color',[0.9 0.9 0.9]);
     hold on;
    end
     
    Fitted=FittedSegment;
    FittedMean=mean(Fitted(:,1));% use this single number when dividing 
    DeltaFoverF=(Raw-Fitted(:,1))/FittedMean; % This is DF/F0 
    

    
    combinedFibVidData(:,6)=DeltaFoverF;% Output Array 
    finalCombinedData(:,4)=DeltaFoverF;
    
    
    sigma=std(DeltaFoverF);
    meanF=mean(DeltaFoverF);
    sigmaplot=(DeltaFoverF-meanF)/sigma;
    
    combinedFibVidData(:,7)=sigmaplot;
    fiberZScores(:,2)=sigmaplot;
    
    
    
    for jj=1:lengthFiberData;
        
        temp=abs(videoDataClean(:,2)-combinedFibVidData(jj,3));
        tempmin=min(temp);
        minindex=find(temp==tempmin);
        vindex=minindex(1);
        
        %edited 11/10/2023
        combinedFibVidData(jj,8)=videoDataClean(vindex,6);
        combinedFibVidData(jj,9)=videoDataClean(vindex,10);
        combinedFibVidData(jj,10)=videoDataClean(vindex,14);
        combinedFibVidData(jj,11)=videoDataClean(vindex,15); 
    end
    
    
    finalCombinedData(:,5:8)=combinedFibVidData(:,8:11);  % interaction zone, corner zone , other 
    
    
  
     
    savefilename0=[basename1 '-ALL.txt'];
    dlmwrite(savefilename0,combinedFibVidData);
    % 
    % savefilename1=[basename1 '-DFF0.txt'];
    % dlmwrite(savefilename1,finalCombinedData);
    % 
    % 
    % savefilename=[basename1 '-sigmaonly.txt'];
    % dlmwrite(savefilename,fiberZScores);
    % 
    MAD=median(abs(DeltaFoverF-meanF));
    [pks1,locs1]=findpeaks(DeltaFoverF-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
    [pks0,locs0]=findpeaks(DeltaFoverF-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);% AUG 2025 PKN
    
    %[pks1,locs1]=findpeaks(sigmaplot,'MinPeakHeight',2,'MinPeakDistance',3);
    numpks=length(locs1);
    numpks0=length(locs0);%AUG 2025 PKN
    
    
    peaktable=zeros(numpks,3);%initialize peaktable 
    peaktable(:,1)=locs1;
    peaktable(:,2)=pks1;
    
    dff0pks1=DeltaFoverF(locs1);
    peaktable(:,3)=dff0pks1;
    peaktable(:,4) = fiberZScores(locs1,2);
    peaktable(:,5) = combinedFibVidData(locs1,3);
    locs2 = combinedFibVidData(locs1,3); % extracts corresponding time for an occurring peak

    %initialize a table containing all signals before threshold %AUG 2025 PKN
    Sigtable=zeros(numpks0,3);
    Sigtable(:,1)=locs0;
    Sigtable(:,2)=pks0;
    dff0pks0=DeltaFoverF(locs0);
    Sigtable(:,3)=dff0pks0;
    Sigtable(:,4) = fiberZScores(locs0,2);
    Sigtable(:,5) = combinedFibVidData(locs0,3);
    locs= combinedFibVidData(locs0,3); % extracts corresponding time for an occurring signal
    
    figure(H0);
    subplot(5,2,2)
    plot(fiberZScores(:,1),fiberZScores(:,2),'Color','black');
    hold on 
    for ii=1:numpks;
        scatter(locs1(ii),pks1(ii),5,'blue');
        hold on
    end
    
    
    figure(H0);
    subplot(5,2,[5,6])
    plot(combinedFibVidData(:,1),combinedFibVidData(:,6),'Color','black');
    hold on 
    for ii=1:numpks;
        scatter(locs1(ii),dff0pks1(ii),5,'blue');
        hold on
    end
    
    counter=1;
    countermax=sum(combinedFibVidData(:,8));
    CCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,8)>0;
        CCC(counter,1)=combinedFibVidData(iii,1);
        CCC(counter,2)=combinedFibVidData(iii,8);
        counter=counter+1;
        else
        end
    end
    plot(CCC(:,1),-0.03*CCC(:,2),'.','Color', [1, 0.156, 0.588]);%pink
    hold on;
    
    counter=1;
    countermax=sum(combinedFibVidData(:,9));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,9)>0;
        CCCC(counter,1)=combinedFibVidData(iii,1);
        CCCC(counter,2)=combinedFibVidData(iii,9);
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.035*CCCC(:,2),'.','Color', [0.145, 0.537, 0.741]);
    ylim([-0.050, max(ylim)]);
    
    
        counter=1;
    countermax=sum(combinedFibVidData(:,10));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,10)>0;
        CCCC(counter,1)=combinedFibVidData(iii,1);
        CCCC(counter,2)=combinedFibVidData(iii,10);
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.040*CCCC(:,2),'.','Color', [0, 1, 0]);
    ylim([-0.050, max(ylim)]);
    
    %second graph with time in s as x-axis (time from FP) %AUG 2025 PKN

    figure(H0);
    subplot(5,2,[7,8])
    xlabel('Time in seconds');

   plot(combinedFibVidData(:,3),combinedFibVidData(:,6),'Color','black');
    hold on 
    for ii=1:numpks;
        scatter(locs2(ii),dff0pks1(ii),5,'blue');
        hold on
    end
    %08/25_PKN


    counter=1;
    countermax=sum(combinedFibVidData(:,8));
    CCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,8)>0;
        CCC(counter,1)=combinedFibVidData(iii,1);
        CCC(counter,2)=combinedFibVidData(iii,3);
        CCC(counter,3)=combinedFibVidData(iii,8);
        counter=counter+1;
        else
        end
    end
    % ✅ Skip plotting if CCC is empty
    if ~isempty(CCC)
    plot(CCC(:,2), -0.03*CCC(:,3), '.', 'Color', [1, 0.156, 0.588]); % pink
    
    end 
    
    %plot(CCC(:,2),-0.03*CCC(:,3),'.','Color', [1, 0.156, 0.588]); %pink
    hold on;
   
    counter=1;
    countermax=sum(combinedFibVidData(:,9));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,9)>0;
        CCCC(counter,1)=combinedFibVidData(iii,1);
        CCCC(counter,2)=combinedFibVidData(iii,3);
        CCCC(counter,3)=combinedFibVidData(iii,9);
        counter=counter+1;
        else
        end
    end
    
    % ✅ Skip plotting if CCC is empty
    if ~isempty(CCCC) 
    plot(CCCC(:,2),-0.035*CCCC(:,3),'.','Color', [0.145, 0.537, 0.741]);
    ylim([-0.050, max(ylim)]);
    end

    counter=1;
    countermax=sum(combinedFibVidData(:,10));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,10)>0;
        CCCC(counter,1)=combinedFibVidData(iii,1);
        CCCC(counter,2)=combinedFibVidData(iii,3);
        CCCC(counter,3)=combinedFibVidData(iii,10);
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,2),-0.040*CCCC(:,3),'.','Color', [0, 1, 0]);
    ylim([-0.050, max(ylim)]);
    
    % ploting a graph showing all signals before threshold with the zone
% 02/28/25PKN
    figure(H0);
    subplot(5,2,[9,10])
    xlabel('Frame Number');
    plot(combinedFibVidData(:,1),combinedFibVidData(:,6),'Color','black');
    hold on 
    for ii=1:numpks0;
        scatter(locs0(ii),dff0pks0(ii),5,'red');
        hold on
    end
    
    counter=1;
    countermax=sum(combinedFibVidData(:,8));
    CCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,8)>0;
        CCC(counter,1)=combinedFibVidData(iii,1);
        CCC(counter,2)=combinedFibVidData(iii,8);
        counter=counter+1;
        else
        end
    end
    plot(CCC(:,1),-0.03*CCC(:,2),'.','Color', [1, 0.156, 0.588]);%pink
    hold on;
    
    counter=1;
    countermax=sum(combinedFibVidData(:,9));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,9)>0;
        CCCC(counter,1)=combinedFibVidData(iii,1);
        CCCC(counter,2)=combinedFibVidData(iii,9);
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.035*CCCC(:,2),'.','Color', [0.145, 0.537, 0.741]);
    ylim([-0.050, max(ylim)]);
   

    counter=1;
    countermax=sum(combinedFibVidData(:,10));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthFiberData;
        if combinedFibVidData(iii,10)>0;
        CCCC(counter,1)=combinedFibVidData(iii,1);
        CCCC(counter,2)=combinedFibVidData(iii,10);
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.040*CCCC(:,2),'.','Color', [0, 1, 0]);
    ylim([-0.050, max(ylim)]);
    
    
    figurefilename2=[basename1 '-.png'];
    saveas(H0,figurefilename2);

    %%
    %creating separate events for averaging

    finalCombinedDataTemp(:,1) = bwlabel(logical(finalCombinedData(:,5)));
    maxNumEvents1 = max(finalCombinedDataTemp(:,1));
    
    if maxNumEvents1 ~=0
        
    TempCellArrayEvent1{maxNumEvents1,13}=[]; 
    
        for i = 1: maxNumEvents1
        locs3 = find(finalCombinedDataTemp(:,1) == i);
        if numel(locs3) > 5;
        FiberDataArray.Active{i} = finalCombinedData(locs3, 4);
        ZFiberDataArray.Active{i} = fiberZScores(locs3, 2);
        finalCombinedDataTemp1 = finalCombinedData(locs3, 4);
        [pksEvent1,locsEvent1]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
        TempCellArrayEvent1{i,1}=basename1;
        TempCellArrayEvent1{i,2}=basename2;
        TempCellArrayEvent1{i,3}="Zone1";
        TempCellArrayEvent1{i,4}=mean(finalCombinedData(locs3, 4));
        TempCellArrayEvent1{i,5}=mean(sigmaplot(locs3));
        TempCellArrayEvent1{i,6}=numel(locsEvent1);
        TempCellArrayEvent1{i,7}=mean(finalCombinedDataTemp1(locsEvent1));
        TempCellArrayEvent1{i,8} = str2num(reply1);
        TempCellArrayEvent1{i,9} =numel(finalCombinedDataTemp1);
        sigmatemp1=sigmaplot(locs3); %02/02/2024 HT
        TempCellArrayEvent1{i,10}=mean(sigmatemp1(locsEvent1)); % 02/02/2024 HT

        [pksEvent01,locsEvent01]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %AUG 2025
        TempCellArrayEvent1{i,11}= numel(locsEvent01);%AUG 2025
        TempCellArrayEvent1{i,12}= mean(finalCombinedDataTemp1(locsEvent01));%AUG 2025
        TempCellArrayEvent1{i,13}=mean(sigmatemp1(locsEvent01));%AUG 2025

        end
        end
    else 
        TempCellArrayEvent1=[];
    end


    finalCombinedDataTemp(:,2) = bwlabel(logical(finalCombinedData(:,6)));
    maxNumEvents2 = max(finalCombinedDataTemp(:,2));
   

    if maxNumEvents2 ~=0
        
     TempCellArrayEvent2{maxNumEvents2,13}=[]; 
     
        for i = 1: maxNumEvents2
        locs3 = find(finalCombinedDataTemp(:,2) == i);
        if numel(locs3) > 5;
        FiberDataArray.Corner{i} = finalCombinedData(locs3, 4);
        ZFiberDataArray.Corner{i} = fiberZScores(locs3, 2);
        finalCombinedDataTemp2 = finalCombinedData(locs3, 4);
        
        [pksEvent2,locsEvent2]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
        TempCellArrayEvent2{i,1}=basename1;
        TempCellArrayEvent2{i,2}=basename2;
        TempCellArrayEvent2{i,3}="Zone2";
        TempCellArrayEvent2{i,4}=mean(finalCombinedData(locs3, 4));
        TempCellArrayEvent2{i,5}=mean(sigmaplot(locs3));
        TempCellArrayEvent2{i,6}=numel(locsEvent2);
        TempCellArrayEvent2{i,7}=mean(finalCombinedDataTemp2(locsEvent2));
        TempCellArrayEvent2{i,8} = str2num(reply1);
        TempCellArrayEvent2{i,9} = numel(finalCombinedDataTemp2);
        sigmatemp2=sigmaplot(locs3); %02/02/2024 HT
        TempCellArrayEvent2{i,10}=mean(sigmatemp2(locsEvent2));

        [pksEvent02,locsEvent02]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%AUG 2025
        TempCellArrayEvent2{i,11}= numel(finalCombinedData(locsEvent02, 4));%AUG 2025
        TempCellArrayEvent2{i,12}= mean(finalCombinedData(locsEvent02, 4));%AUG 2025
        TempCellArrayEvent2{i,13}=mean(sigmatemp2(locsEvent02));%AUG 2025
        end
        end
    else 
    TempCellArrayEvent2=[];     
    end

    
    finalCombinedDataTemp(:,3) = bwlabel(logical(finalCombinedData(:,7)));
    maxNumEvents3 = max(finalCombinedDataTemp(:,3));
    

    if maxNumEvents3 ~=0
        
    TempCellArrayEvent3{maxNumEvents3,13}=[]; 
    
        for i = 1: maxNumEvents3
        locs3 = find(finalCombinedDataTemp(:,3) == i);
        FiberDataArray.NonActive{i} = finalCombinedData(locs3, 4);
        ZFiberDataArray.NonActive{i} = fiberZScores(locs3, 2);
        if numel(locs3) > 5
            finalCombinedDataTemp3 = finalCombinedData(locs3, 4);
            [pksEvent3,locsEvent3]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
            TempCellArrayEvent3{i,1}=basename1;
            TempCellArrayEvent3{i,2}=basename2;
            TempCellArrayEvent3{i,3}="Zone3";
            TempCellArrayEvent3{i,4}=mean(finalCombinedData(locs3, 4));
            TempCellArrayEvent3{i,5}=mean(sigmaplot(locs3));
            TempCellArrayEvent3{i,6}=numel(locsEvent3);
            TempCellArrayEvent3{i,7}=mean(finalCombinedDataTemp3(locsEvent3));
            TempCellArrayEvent3{i,8} = str2num(reply1);
            TempCellArrayEvent3{i,9} =numel(finalCombinedDataTemp3);
            sigmatemp3=sigmaplot(locs3); %02/02/2024 HT
            TempCellArrayEvent3{i,10}=mean(sigmatemp3(locsEvent3));
            
            [pksEvent03,locsEvent03]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%AUG 2025
            TempCellArrayEvent3{i,11}= numel(finalCombinedData(locsEvent03, 4));%AUG 2025
            TempCellArrayEvent3{i,12}= mean(finalCombinedData(locsEvent03, 4));%AUG 2025
            TempCellArrayEvent3{i,13}=mean(sigmatemp3(locsEvent03));%AUG 2025
        end
        end
    else 
        
    TempCellArrayEvent3=[];    
    end


    
    finalCombinedDataTemp(:,4) = bwlabel(logical(finalCombinedData(:,8)));
    maxNumEvents4 = max(finalCombinedDataTemp(:,4));
    

    if maxNumEvents4 ~=0
        
    TempCellArrayEvent4{maxNumEvents4,13}=[]; 
    
        for i = 1: maxNumEvents4
        locs4 = find(finalCombinedDataTemp(:,4) == i);
        FiberDataArray.NonZone{i} = finalCombinedData(locs4, 4);
        ZFiberDataArray.NonZone{i} = fiberZScores(locs4, 2);
        if numel(locs4) > 5
            finalCombinedDataTemp4 = finalCombinedData(locs4, 4);
            [pksEvent4,locsEvent4]=findpeaks(finalCombinedData(locs4, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
            TempCellArrayEvent4{i,1}=basename1;
            TempCellArrayEvent4{i,2}=basename2;
            TempCellArrayEvent4{i,3}="OtherZone";
            TempCellArrayEvent4{i,4}=mean(finalCombinedData(locs4, 4));
            TempCellArrayEvent4{i,5}=mean(sigmaplot(locs4));
            TempCellArrayEvent4{i,6}=numel(locsEvent4);
            TempCellArrayEvent4{i,7}=mean(finalCombinedDataTemp4(locsEvent4));
            TempCellArrayEvent4{i,8} = str2num(reply1);
            TempCellArrayEvent4{i,9} =numel(finalCombinedDataTemp4);
            sigmatemp4=sigmaplot(locs4); %02/02/2024 HT
            TempCellArrayEvent4{i,10}=mean(sigmatemp4(locsEvent4));%02/02/2024 HT

            [pksEvent04,locsEvent04]=findpeaks(finalCombinedData(locs4, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%AUG 2025
            TempCellArrayEvent4{i,11}= numel(finalCombinedData(locsEvent04, 4));%AUG 2025
            TempCellArrayEvent4{i,12}= mean(finalCombinedData(locsEvent04, 4));%AUG 2025
            TempCellArrayEvent4{i,13}=mean(sigmatemp4(locsEvent04));%AUG 2025
        end
        end
    else 
        
    TempCellArrayEvent4=[];    
    end
    
    
    
    
   TempCellArrayAllEvents = [TempCellArrayEvent1 ; TempCellArrayEvent2; TempCellArrayEvent3; TempCellArrayEvent4]; %modified 11/10/23

   TempCellLabelEvents={"basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z", "numSig" "meanSigDFF" "MeanSigDFF in Z"};
   TempCellLabelEvents1=["basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z", "numSig" "meanSigDFF" "MeanSigDFF in Z"];


   TableTempCellLabelEvents=cell2table(TempCellLabelEvents, "VariableNames", TempCellLabelEvents1);
   datafilenameEvents=['125_d-SLR_SAMPLE_events_all_081325_PKN' '.csv'];
    if isfile(datafilenameEvents);
      
    else
    writetable(TableTempCellLabelEvents,datafilenameEvents);
    end
    
     CellTableEvents=table2cell(readtable(datafilenameEvents));
     CellTableNewEvents=[CellTableEvents;TempCellArrayAllEvents];
     TableNewEvents=cell2table(CellTableNewEvents, "VariableNames", TempCellLabelEvents1);
     writetable(TableNewEvents,datafilenameEvents);


    datafilenameEventsMat = [basename1 'events.mat'];

    %matlab file for all 3 events in a data structure
    save(datafilenameEventsMat, "FiberDataArray");

    
    %calculate average fluorescence during the zone 
    
    TotalZone1Count=sum(finalCombinedData(:,5));
    TempArray=zeros(TotalZone1Count,1);
    
    TotalZone2Count=sum(finalCombinedData(:,6));
    TempArray2=zeros(TotalZone2Count,1);
    
    TotalZone3Count=sum(finalCombinedData(:,7));
    TempArray3=zeros(TotalZone3Count,1);
    
    TotalOtherZoneCount=sum(finalCombinedData(:,8));
    TempArray4=zeros(TotalOtherZoneCount,1);
    


    
    nn=1; %counter 
    nnn=1; %counter 
    nnnn=1; %counter 
    nnnnn=1; % counter 11/10/23
    
    for iii=1:lengthFiberData 
     if finalCombinedData(iii,5)==1;
         TempArray(nn,1)=finalCombinedData(iii,4);
         nn=nn+1;
     elseif finalCombinedData(iii,6)==1
         TempArray2(nnn,1)=finalCombinedData(iii,4);
         nnn=nnn+1;
     elseif finalCombinedData(iii,7)==1
         TempArray3(nnnn,1)=finalCombinedData(iii,4);
         nnnn=nnnn+1;
     else
         TempArray4(nnnnn,1)=finalCombinedData(iii,4);
         nnnnn=nnnnn+1;
     end
    end

    %Finding DF/F signal amplitude (before thresold) %AUG 2025 PKN 
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks0;
        NNN=Sigtable(iii,1);
        if finalCombinedData(NNN,5)==1;
            Zone1Sig(nn,1)=Sigtable(iii,3);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            Zone2Sig(nnn,1)=Sigtable(iii,3);
            nnn=nnn+1;
        elseif finalCombinedData(NNN,7)==1;
            Zone3Sig(nnnn,1)=Sigtable(iii,3);
            nnnn=nnnn+1;
        else
            OtherZoneSig(nnnnn,1)=Sigtable(iii,3);
            nnnnn=nnnnn+1;
        end 
    end
        
  %Finding Z-score of DF/F signal amplitude (before thresold) %AUG2025 PKN 
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks0;
        NNN=Sigtable(iii,1);
        if finalCombinedData(NNN,5)==1;
            ZscoreZone1Sig(nn,1)=Sigtable(iii,4);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            ZscoreZone2Sig(nnn,1)=Sigtable(iii,4);
            nnn=nnn+1;
        elseif finalCombinedData(NNN,7)==1
            ZscoreZone3Sig(nnnn,1)=Sigtable(iii,4);
            nnnn=nnnn+1;
        else
            ZscoreOtherZoneSig(nnnnn,1)=Sigtable(iii,4);
            nnnnn=nnnnn+1;
        end 
    end


    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    nnnnn=1; 
    for iii=1:numpks;
        NNN=peaktable(iii,1);
        if finalCombinedData(NNN,5)==1;
            Zone1Peak(nn,1)=peaktable(iii,3);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            Zone2Peak(nnn,1)=peaktable(iii,3);
            nnn=nnn+1;
        elseif finalCombinedData(NNN,7)==1;
            Zone3Peak(nnnn,1)=peaktable(iii,3);
            nnnn=nnnn+1;
        else
            OtherZonePeak(nnnnn,1)=peaktable(iii,3);
            nnnnn=nnnnn+1;
        end 
    end
   
    % Finding Z-score of peak amplitude  %AUG 2025 PKN     
  nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;

    for iii=1:numpks;
        NNN=peaktable(iii,1);
        if finalCombinedData(NNN,5)==1;
            ZscoreZone1Peak(nn,1)=peaktable(iii,4);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            ZscoreZone2Peak(nnn,1)=peaktable(iii,4);
            nnn=nnn+1;
        elseif finalCombinedData(NNN,7)==1
            ZscoreZone3Peak(nnnn,1)=peaktable(iii,4);
            nnnn=nnnn+1;
        else
            ZscoreOtherZonePeak(nnnnn,1)=peaktable(iii,3);
            nnnnn=nnnnn+1;
        end 
    end  
    
   if nn==1;
       Zone1Peak=[];
   end
   
   if nnn==1;
       Zone2Peak=[];
   else
   end
   
   if nnnn==1;
       Zone3Peak=[];
   else
   end
   
   if nnnnn==1;
       OtherZonePeak=[];
   else
   end
   if nn==1;
       ZscoreZone1Peak=[];
   end
   
   if nnn==1;
       ZscoreZone2Peak=[];
   else
   end
   
   if nnnn==1
       ZscoreZone3Peak=[];
   else
   end
   if nnnnn==1;
       ZscoreOtherZonePeak=[];
   else
   end

   Zone1Sig=[]
    if ~isempty(Zone1Sig)
    % The array has data, so calculate the mean.
      DFFZone1SigMean= mean(Zone1Sig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
      DFFZone1SigMean = NaN; 
    end

    ZscoreZone1Sig=[]
    if ~isempty(ZscoreZone1Sig)
    % The array has data, so calculate the mean.
      ZscoreDFFZone1SigMean= mean(ZscoreZone1Sig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
      ZscoreDFFZone1SigMean = NaN; 
    end

   DFFZone1Mean=mean(TempArray);
   DFFZone1SigMean=mean(Zone1Sig);
   ZscoreDFFZone1SigMean=mean(ZscoreZone1Sig);
   DFFZone1Stdev=std(TempArray);
   Zone1Counts=size(TempArray,1);
   DFFZone1SigCounts=size(Zone1Sig,1);
   Zone1SignalRate=DFFZone1SigCounts/Zone1Counts*20;% peak#/sec

 
   DFFZone2Mean=mean(TempArray2);
   DFFZone2SigMean=mean(Zone2Sig);
   ZscoreDFFZone2SigMean=mean(ZscoreZone2Sig);
   DFFZone2Stdev=std(TempArray2);
   Zone2Counts=size(TempArray2,1);
   DFFZone2SigCounts=size(Zone2Sig,1);
   Zone2SignalRate=DFFZone2SigCounts/Zone2Counts*20;% peak#/sec
    

   DFFZone3Mean=mean(TempArray3);
   DFFZone3SigMean=mean(Zone3Sig);
   ZscoreDFFZone3SigMean=mean(ZscoreZone3Sig);
   DFFZone3Stdev=std(TempArray3);
   Zone3Counts=size(TempArray3,1); 
   DFFZone3SigCounts=size(Zone3Sig,1);
   Zone3SignalRate=DFFZone3SigCounts/Zone3Counts*20;% peak#/sec
   
   
   DFFOtherZoneMean=mean(TempArray4);
   DFFOtherZoneSigMean=mean(OtherZoneSig);
   ZscoreDFFOtherZoneSigMean=mean(ZscoreOtherZoneSig);
   DFFOtherZoneStdev=std(TempArray4);
   OtherZoneCounts=size(TempArray4,1); 
   DFFOtherZoneSigCounts=size(OtherZoneSig,1);
   OtherZoneSignalRate=DFFOtherZoneSigCounts/OtherZoneCounts*20;% peak#/sec
    
   Zone1PeakCounts=size(Zone1Peak,1);
   RateZone1PeakCounts=Zone1PeakCounts/Zone1Counts*20;% peak#/sec
   Zone1PeakHeight=mean(Zone1Peak);
   ZscoreZone1PeakHeight=mean(ZscoreZone1Peak);
    
   Zone2PeakCounts=size(Zone2Peak,1);
   RateZone2PeakCounts=Zone2PeakCounts/Zone2Counts*20;% peak#/sec
   Zone2PeakHeight=mean(Zone2Peak);
   ZscoreZone2PeakHeight=mean(ZscoreZone2Peak);
    
   Zone3PeakCounts=size(Zone3Peak,1);
   RateZone3PeakCounts=Zone3PeakCounts/Zone3Counts*20;% peak#/sec
   Zone3PeakHeight=mean(Zone3Peak);
   ZscoreZone3PeakHeight=mean(ZscoreZone3Peak);
   
   OtherZonePeakCounts=size(OtherZonePeak,1);
   RateOtherZonePeakCounts=OtherZonePeakCounts/OtherZoneCounts*20;% peak#/sec
   OtherZonePeakHeight=mean(OtherZonePeak); 
   ZscoreOtherZonePeakHeight=mean(ZscoreOtherZonePeak);
   
   
   TempCellArray{1,50}=[]; 
   TempCellArray{1,1}=basename1;
   TempCellArray{1,2}=basename2;
   
   TempCellArray{1,3}=NSegment;
   
   TempCellArray{1,4}=startpoint;
   TempCellArray{1,5}=endpoint;
  
   TempCellArray{1,6}=DFFZone1Mean;
   TempCellArray{1,7}=DFFZone1SigMean;
   TempCellArray{1,8}=ZscoreDFFZone1SigMean;
   TempCellArray{1,9}=DFFZone1SigCounts;
   TempCellArray{1,10}=Zone1SignalRate;
   TempCellArray{1,11}=DFFZone1Stdev;
   TempCellArray{1,12}=Zone1Counts;
   
   TempCellArray{1,13}=DFFZone2Mean;
   TempCellArray{1,14}=DFFZone2SigMean;
   TempCellArray{1,15}=ZscoreDFFZone2SigMean;
   TempCellArray{1,16}=DFFZone2SigCounts;
   TempCellArray{1,17}=Zone2SignalRate;
   TempCellArray{1,18}=DFFZone2Stdev;
   TempCellArray{1,19}=Zone2Counts;
   
   TempCellArray{1,20}=DFFZone3Mean;
   TempCellArray{1,21}=DFFZone3SigMean;
   TempCellArray{1,22}=ZscoreDFFZone3SigMean;
   TempCellArray{1,23}=DFFZone3SigCounts;
   TempCellArray{1,24}=Zone3SignalRate;
   TempCellArray{1,25}=DFFZone3Stdev;
   TempCellArray{1,26}=Zone3Counts;

   TempCellArray{1,27}=DFFOtherZoneMean;
   TempCellArray{1,28}=DFFOtherZoneSigMean;
   TempCellArray{1,29}=ZscoreDFFOtherZoneSigMean;
   TempCellArray{1,30}=DFFOtherZoneSigCounts;
   TempCellArray{1,31}=OtherZoneSignalRate;
   TempCellArray{1,32}=DFFOtherZoneStdev;
   TempCellArray{1,33}=OtherZoneCounts;

   TempCellArray{1,34}=Zone1PeakHeight;
   TempCellArray{1,35}=ZscoreZone1PeakHeight;
   TempCellArray{1,36}=RateZone1PeakCounts;
   TempCellArray{1,37}=Zone1PeakCounts;

   TempCellArray{1,38}=Zone2PeakHeight;
   TempCellArray{1,39}=ZscoreZone2PeakHeight;
   TempCellArray{1,40}=RateZone2PeakCounts;
   TempCellArray{1,41}=Zone2PeakCounts;

   TempCellArray{1,42}=Zone3PeakHeight;
   TempCellArray{1,43}=ZscoreZone3PeakHeight;
   TempCellArray{1,44}=RateZone3PeakCounts;
   TempCellArray{1,45}=Zone3PeakCounts;

   TempCellArray{1,46}=OtherZonePeakHeight;
   TempCellArray{1,47}=ZscoreOtherZonePeakHeight;
   TempCellArray{1,48}=RateOtherZonePeakCounts;
   TempCellArray{1,49}=OtherZonePeakCounts;
      
   
   TempCellArray{1,50}=reply1;
 
    
%    TempCellLabel={'basename1' 'basename' 'NSegment' 'datastart' 'dataend' 'DFFZoneMean' 'DFFZoneStdev' 'DFFZoneCounts' 'DFFCornerZoneMean' 'DFFCornerZoneStdev' 'DFFCornerZoneCounts' 'DFFNonZoneMean' 'DFFNonZoneStdev' 'DFFNonZoneCounts' 'ZonePeakCounts' 'NormalizedZonePeakCounts' 'ZonePeakHeight' 'CornerZonePeakCounts' 'NormalizedCornerZonePeakCounts' 'CornerZonePeakHeight' 'NonZonePeakCounts' 'NormalizedNonZonePeakCounts' 'NonZonePeakHeight'};
   
% TempCellLabel={"basename1" "basename" "NSegment" "datastart" "dataend" "DFFZone1Mean" "DFFZone1Stdev" "DFFZone1Counts" "DFFCornerZoneMean" "DFFCornerZoneStdev" "DFFCornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneStdev" "DFFNonZoneCounts" "ZonePeakCounts" "NormalizedZonePeakCounts" "ZonePeakHeight" "CornerZonePeakCounts" "NormalizedCornerZonePeakCounts" "CornerZonePeakHeight" "NonZonePeakCounts" "NormalizedNonZonePeakCounts" "NonZonePeakHeight" "TypeOfAnalysisControl"};

TempCellLabel={"basename1" "basename2" "NSegment" "startpoint" "endpoint" "DFFZone1Mean" "DFFZone1SigMean" "ZscoreDFFZone1SigMean" "DFFZone1SigCounts" "Zone1SignalRate" "DFFZone1Stdev" "Zone1Counts" "DFFZone2Mean" "DFFZone2SigMean" "ZscoreDFFZone2SigMean" "DFFZone2SigCounts" "Zone2SignalRate" "DFFZone2Stdev" "Zone2Counts" "DFFZone3Mean" "DFFZone3SigMean" "ZscoreDFFZone3SigMean" "DFFZone3SigCounts" "Zone3SignalRate" "DFFZone3Stdev" "Zone3Counts" "DFFOtherZoneMean" "DFFOtherZoneSigMean" "ZscoreDFFOtherZoneSigMean" "DFFOtherZoneSigCounts" "OtherZoneSignalRate" "DFFOtherZoneStdev" "OtherZoneCounts" "Zone1PeakHeight" "ZscoreZone1PeakHeight" "RateZone1PeakCounts" "Zone1PeakCounts" "Zone2PeakHeight"  "ZscoreZone2PeakHeight" "RateZone2PeakCounts" "Zone2PeakCounts" "Zone3PeakHeight" "ZscoreZone3PeakHeight" "RateZone3PeakCounts" "Zone3PeakCounts"  "OtherZonePeakHeight" "ZscoreOtherZonePeakHeight" "RateOtherZonePeakCounts" "OtherZonePeakCounts" "Type of Analysis Control"} 
            

TempCellLabel1=["basename1" "basename2" "NSegment" "startpoint" "endpoint" "DFFZone1Mean" "DFFZone1SigMean" "ZscoreDFFZone1SigMean" "DFFZone1SigCounts" "Zone1SignalRate" "DFFZone1Stdev" "Zone1Counts" "DFFZone2Mean" "DFFZone2SigMean" "ZscoreDFFZone2SigMean" "DFFZone2SigCounts" "Zone2SignalRate" "DFFZone2Stdev" "Zone2Counts" "DFFZone3Mean" "DFFZone3SigMean" "ZscoreDFFZone3SigMean" "DFFZone3SigCounts" "Zone3SignalRate" "DFFZone3Stdev" "Zone3Counts" "DFFOtherZoneMean" "DFFOtherZoneSigMean" "ZscoreDFFOtherZoneSigMean" "DFFOtherZoneSigCounts" "OtherZoneSignalRate" "DFFOtherZoneStdev" "OtherZoneCounts" "Zone1PeakHeight" "ZscoreZone1PeakHeight" "RateZone1PeakCounts" "Zone1PeakCounts" "Zone2PeakHeight"  "ZscoreZone2PeakHeight" "RateZone2PeakCounts" "Zone2PeakCounts" "Zone3PeakHeight" "ZscoreZone3PeakHeight" "RateZone3PeakCounts" "Zone3PeakCounts"  "OtherZonePeakHeight" "ZscoreOtherZonePeakHeight" "RateOtherZonePeakCounts" "OtherZonePeakCounts" "Type of Analysis Control"] 

%    datafilename=['FiberDataExport4' '.csv'];
%     if ~isempty(datafilename);
%         xlswrite(datafilename,TempCellLabel,'Sheet1');
%     else
%     end
% 
%     [success,message]=xlsappend(datafilename,TempCellArray);
   TableTempCellLabel=cell2table(TempCellLabel, "VariableNames", TempCellLabel1);
   datafilename=['125_d-SLR_SAMPLE_summary_all_081325_PKN' '.csv'];
    if isfile(datafilename);
      
    else
    writetable(TableTempCellLabel,datafilename);
    end
    
     CellTable=table2cell(readtable(datafilename));
     CellTableNew=[CellTable;TempCellArray];
     TableNew=cell2table(CellTableNew, "VariableNames", TempCellLabel1);
     writetable(TableNew,datafilename);




%%
   reply = input('Do you want partial analysis? y/n [y]:','s');
   
   if isempty(reply);
      reply = 'n';
   end
   
   if reply =='y';
       startpoint = input('starting datapoint  ');
       endpoint = input('endding datapoint  ');
%    else
%    end
   startpoint=int16(startpoint); 
   endpoint=int16(endpoint); 
   
   partialAnalysisData=[];
   partialAnalysisData=zeros((endpoint-startpoint+1),9);
   partialAnalysisData(:,1:3)=finalCombinedData(startpoint:endpoint,1:3);
   partialAnalysisData(:,4:5)=fiberData(startpoint:endpoint,4:5);
   partialAnalysisData(:,7:10)=finalCombinedData(startpoint:endpoint,5:8); %11/17/2023 8 column - 6??
   
   lengthF=size(partialAnalysisData,1);

    tempControlPA=partialAnalysisData(:,4);

    RawPA=partialAnalysisData(:,5);

   if isempty(reply1)
      reply1 = '1';
   end
   
   if reply1 =='1'
       ControlPA = smoothdata(tempControlPA, 'rlowess', 5);
   elseif reply1 == '2'
       ControlPA = smoothdata(tempControlPA, 'rlowess', 500);
   elseif reply1 == '3'
%        ControlPA = smoothdata(RawPA, 'rlowess', 500);
 ControlPA = sgolayfilt(RawPA, 3 , 501);
   end
    
    
    SegmentLengthF=floor(lengthF/NSegment);
    SubtractedF=zeros(lengthF,1);%% Preallocate
    FittedSegmentF=zeros(lengthF,1);%% Preallocate
%     SegmentFrames=zeros(NSegment,1);
    
    SegmentPoints=ones(NSegment,1);

    for ii=1:NSegment;
    data415=[];
    data470=[];
    
        if ii ~= NSegment;
    mm = 1+SegmentLengthF*(ii-1);
    nn = SegmentLengthF*(ii);
    SegmentPoints(ii,1)=mm;
         else
    mm = 1+SegmentLengthF*(ii-1);
    nn= lengthF;
    SegmentPoints(ii,1)=mm;
        end
        
    data415=ControlPA(mm:nn,1);
    data470=RawPA(mm:nn,1);
    

%     IniValues=[1.000 -0.002];
%     options = optimset('Display','notify','MaxIter',1000000,'MaxFunEvals',1000000);
%      Estimates = fminsearch(@myfitFPSY,IniValues,options,TempC, TempA);
%     Estimates = fminsearchbnd(@myfitFPSY,IniValues,[0.1000 -inf],[5.000 inf],options,TempC, TempA);
%     
%     %%%Results 
%     a=Estimates(1);
%     b=Estimates(2);

    P = polyfit(data415, data470, 3);

    %%%Subtract Iso from Calcium 
    TempFittedF=P(1)*data415.^3 + P(2)*data415.^2 + P(3)*data415 + P(4);
    TempSubtractedF=(data470 - TempFittedF);

    %%%Subtract Iso from Calcium    
    SubtractedA(mm:nn)=TempSubtractedF;
    FittedSegmentF(mm:nn)=TempFittedF;
  
            
    end
    
    
    H1=figure('Name',[basename1 '_' num2str(startpoint) '_' num2str(endpoint)]);
  
    figure(H1);
    subplot(3,2,1);
    plot(combinedFibVidData(:,4),'Color','blue');
    hold on
    plot(combinedFibVidData(:,5),'Color','green');
    hold on
    plot(Control,'Color','magenta');
    hold on
    minh1=min(min(combinedFibVidData(:,4)), min(combinedFibVidData(:,5)));
    maxh1=max(max(combinedFibVidData(:,4)), max(combinedFibVidData(:,5)));
    xxh1=[startpoint startpoint];
    yyh1=[minh1 maxh1];
    xxh1_=[endpoint endpoint];
    yyh1_=[minh1 maxh1];
    plot(xxh1,yyh1,'Color',[0.7 0.7 0.7]);
    hold on
    plot(xxh1_,yyh1_,'Color',[0.7 0.7 0.7]);
    
    subplot(3,2,[3,4])
    plot(FittedSegmentF,'LineWidth', 2,'Color','green');
    hold on;
    plot(partialAnalysisData(:,5),'Color',[0.4 0.4 0.4]);
    hold on
    for kk=1:NSegment;
     XX=[SegmentPoints(kk,1)  SegmentPoints(kk,1)];
     YY=[min(Raw)  max(Raw)];
     plot(XX,YY,'Color',[0.9 0.9 0.9]);
     hold on;
    end
     
    FittedF=FittedSegmentF;
    FittedMeanF=mean(FittedF(:,1));% use this single number when dividing 
    DeltaFoverFF=(partialAnalysisData(:,5)-FittedF(:,1))/FittedMeanF; % This is DF/F0 
    
%     figure(H1);
%     subplot(3,2,2)
%     plot(DeltaFoverFF);
    
    
    %
    partialAnalysisData(:,6)=DeltaFoverFF;
   
   
    sigmaF=std(DeltaFoverFF);
    meanFF=mean(DeltaFoverFF);
    sigmaplotF=(DeltaFoverFF-meanFF)/sigmaF;
    
    figure(H1);
    subplot(3,2,2)
    plot(sigmaplotF);
    
    
    savefilename0=[basename1 '_' num2str(startpoint) '_' num2str(endpoint) '-ALL.txt'];
    dlmwrite(savefilename0,partialAnalysisData);
 
    MAD=median(abs(DeltaFoverFF-meanFF));
    [pks2,locs2]=findpeaks(DeltaFoverFF-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
    [pks02,locs02]=findpeaks(DeltaFoverFF-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%defining "no thresold" % AUG 2025 PKN
    
    %[pks2,locs2]=findpeaks(sigmaplotF,'MinPeakHeight',2,'MinPeakDistance',3);
    numpks2=length(locs2);
    numpks02=length(locs02);% AUG 2025 PKN
    
       
    peaktable2=zeros(numpks2,3);%initialize peaktable 
    peaktable2(:,1)=locs2;
    peaktable2(:,2)=pks2;
    
    dff0pks2=DeltaFoverFF(locs2);
    peaktable2(:,3)=dff0pks2;
    peaktable2(:,4) = fiberZScores(locs2,2);% AUG 2025 PKN

    %initialize a table containing all signals before threshold % AUG 2025 PKN
    PASigtable=zeros(numpks02,3);
    PASigtable(:,1)=locs02;
    PASigtable(:,2)=pks02;
    dff0pks02=DeltaFoverF(locs02);
    PASigtable(:,3)=dff0pks02;
    PASigtable(:,4) = fiberZScores(locs02,2);
    PASigtable(:,5) = combinedFibVidData(locs02,3);
    locs= combinedFibVidData(locs0,3); % extracts corresponding time for an occurring signal
    
    
    figure(H1);
    subplot(3,2,[5,6])
    plot(partialAnalysisData(:,6),'Color','black');
    hold on 
    for ii=1:numpks2;
        scatter(locs2(ii),dff0pks2(ii),5,'blue');
        hold on
    end
    
    counter=1;
    countermax=sum(partialAnalysisData(:,7));
    CCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthF;
        if partialAnalysisData(iii,7)>0;
        CCC(counter,1)=iii;
        CCC(counter,2)=partialAnalysisData(iii,7);
        counter=counter+1;
        else
        end
    end
    plot(CCC(:,1),-0.03*CCC(:,2),'.','Color',[1, 0.156, 0.588]);
    hold on;
    
    
    
    counter=1;
    countermax=sum(partialAnalysisData(:,8));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthF;
        if partialAnalysisData(iii,8)>0;
        CCCC(counter,1)=iii;
        CCCC(counter,2)=partialAnalysisData(iii,8);
        counter=counter+1;31
        else
        end
    end
    plot(CCCC(:,1),-0.035*CCCC(:,2),'.','Color',[0.145, 0.537, 0.741]);
    
    hold on;
    
    
    counter=1;
    countermax=sum(partialAnalysisData(:,9));
    CCCC=zeros(countermax,2);%subarray for zone only
    for iii=1:lengthF;
        if partialAnalysisData(iii,9)>0;
        CCCC(counter,1)=iii;
        CCCC(counter,2)=partialAnalysisData(iii,9);
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.040*CCCC(:,2),'.','Color',[0, 1, 0]);
    
    
    
    
    
    
    
    
    ylim([-0.050, max(ylim)+0.005]);
    
    
    figurefilename2=[basename1 '_' num2str(startpoint) '_' num2str(endpoint) '-.png'];
    saveas(H1,figurefilename2);
    
    finalCombinedDataTempPA = [];

    
    
    locs3=[];
    %creating separate events for averaging (CHANGE ALL OF THE VARIABLES)

   finalCombinedDataTempPA(:,1) = bwlabel(logical(partialAnalysisData(:,7)));
   maxNumEvents1 = max(finalCombinedDataTempPA(:,1));
   
    if maxNumEvents1 > 0
    TempCellArrayEventPA1=[];    
    TempCellArrayEventPA1{maxNumEvents1,13}=[];     
        for i = 1: maxNumEvents1
        locs3 = find(finalCombinedDataTempPA(:,1) == i);
        if numel(locs3) > 5;
        FiberDataArrayPA.Active{i} = partialAnalysisData(locs3, 4);
        finalCombinedDataTempPA1 = partialAnalysisData(locs3, 4);
        [pksPAEvent1,locsPAEvent1]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);%May23,2023
        TempCellArrayEventPA1{i,1}=basename1;
        TempCellArrayEventPA1{i,2}=basename2;
        TempCellArrayEventPA1{i,3}="Zone1";
        TempCellArrayEventPA1{i,4}=mean(partialAnalysisData(locs3, 6));%May23,2023
        TempCellArrayEventPA1{i,5}=mean(sigmaplotF(locs3)); % from line 802  % mean Z-score 
        TempCellArrayEventPA1{i,6}=numel(locsPAEvent1);
        TempCellArrayEventPA1{i,7}=mean(finalCombinedDataTempPA1(locsPAEvent1));
        TempCellArrayEventPA1{i,8} = str2num(reply1);
        TempCellArrayEventPA1{i,9}=numel(finalCombinedDataTempPA1);
        sigmaplotFtemp1=sigmaplotF(locs3); %HT 02/02/2024
        TempCellArrayEventPA1{i,10}=mean(sigmaplotFtemp1(locsPAEvent1)); % Mean Peak Height Z-Score %02/02/2024

        [pksPAEvent01,locsPAEvent01]=findpeaks(partialAnalysisData(locs3,6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%02/27/2025 PKN
        TempCellArrayEventPA1{i,11}= numel(locsPAEvent01);%02/27/2025 PKN
        TempCellArrayEventPA1{i,12}= mean(finalCombinedDataTempPA1(locsPAEvent01));%03/03/2025 PKN
        TempCellArrayEventPA1{i,13}=mean(sigmaplotFtemp1(locsPAEvent01));%03/03/2025 PKN

                
        end
        end
    else 
        
        TempCellArrayEventPA1=[]; 
    end


    finalCombinedDataTempPA(:,2) = bwlabel(logical(partialAnalysisData(:,8)));
    maxNumEvents2 = max(finalCombinedDataTempPA(:,2));
    
    if maxNumEvents2 > 0
    TempCellArrayEventPA2=[];       
    TempCellArrayEventPA2{maxNumEvents2,10}=[];     
        for i = 1: maxNumEvents2
        locs3 = find(finalCombinedDataTempPA(:,2) == i);
        if numel(locs3) > 5;
        % need to fix the location index to use find
        FiberDataArrayPA.Corner{i} = partialAnalysisData(locs3, 4);
        finalCombinedDataTempPA2 = partialAnalysisData(locs3, 4);
        [pksPAEvent2,locsPAEvent2]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);%May23,2023
        TempCellArrayEventPA2{i,1}=basename1;
        TempCellArrayEventPA2{i,2}=basename2;
        TempCellArrayEventPA2{i,3}="Zone2";
        TempCellArrayEventPA2{i,4}=mean(partialAnalysisData(locs3, 6)); %May 2023
        TempCellArrayEventPA2{i,5}=mean(sigmaplotF(locs3)); % from line 802  % mean Z-score
        TempCellArrayEventPA2{i,6}=numel(locsPAEvent2);
        TempCellArrayEventPA2{i,7}=mean(finalCombinedDataTempPA2(locsPAEvent2));
        TempCellArrayEventPA2{i,8} = str2num(reply1);
        TempCellArrayEventPA2{i,9}=numel(finalCombinedDataTempPA2);
        sigmaplotFtemp2=sigmaplotF(locs3); %HT 02/02/2024
        TempCellArrayEventPA2{i,10}=mean(sigmaplotFtemp2(locsPAEvent2)); % Mean Peak Height Z-Score

        [pksPAEvent02,locsPAEvent02]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%02/27/2025 PKN
        TempCellArrayEventPA2{i,11}= numel(locsPAEvent02);%02/27/2025 PKN
        TempCellArrayEventPA2{i,12}= mean(finalCombinedDataTempPA2(locsPAEvent02));%03/03/2025 PKN
        TempCellArrayEventPA2{i,13}=mean(sigmaplotFtemp2(locsPAEvent02));%03/03/2025 PKN


        end
        end
    else 
        
    TempCellArrayEventPA2=[];     
    end

    
      
    
    finalCombinedDataTempPA(:,3) = bwlabel(logical(partialAnalysisData(:,9)));
    maxNumEvents3 = max(finalCombinedDataTempPA(:,3));
    

    if maxNumEvents3 > 0
    TempCellArrayEventPA3=[];      
    TempCellArrayEventPA3{maxNumEvents3,13}=[];     
        for i = 1: maxNumEvents3
        locs3 = find(finalCombinedDataTempPA(:,3) == i);
         if numel(locs3) > 5
            FiberDataArrayPA.NonActive{i} = partialAnalysisData(locs3, 4);
            finalCombinedDataTempPA3 = partialAnalysisData(locs3, 4);
            [pksPAEvent3,locsPAEvent3]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);%May23,2023
            TempCellArrayEventPA3{i,1}=basename1;
            TempCellArrayEventPA3{i,2}=basename2;
            TempCellArrayEventPA3{i,3}="Zone3";
            TempCellArrayEventPA3{i,4}=mean(partialAnalysisData(locs3, 6));%May23,2023
            TempCellArrayEventPA3{i,5}=mean(sigmaplotF(locs3));% from line 802  % mean Z-score 
            TempCellArrayEventPA3{i,6}=numel(locsPAEvent3);
            TempCellArrayEventPA3{i,7}=mean(finalCombinedDataTempPA3(locsPAEvent3));
            TempCellArrayEventPA3{i,8} =str2num(reply1);
            TempCellArrayEventPA3{i,9}=numel(finalCombinedDataTempPA3);
            sigmaplotFtemp3=sigmaplotF(locs3); %HT 02/02/2024
            TempCellArrayEventPA3{i,10}=mean(sigmaplotFtemp3(locsPAEvent3)); % Mean Peak Height Z-Score

            [pksPAEvent03,locsPAEvent03]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %02/27/2025 PKN
            TempCellArrayEventPA3{i,11}= numel(locsPAEvent03);%02/27/2025 PKN
            TempCellArrayEventPA3{i,12}= mean(finalCombinedDataTempPA3(locsPAEvent03));%03/03/2025 PKN
            TempCellArrayEventPA3{i,13}=mean(sigmaplotFtemp3(locsPAEvent03));%03/03/2025 PKN

         end
        end
    else 
    TempCellArrayEventPA3=[];     
    end

    
    finalCombinedDataTempPA(:,4) = bwlabel(logical(partialAnalysisData(:,10)));
    maxNumEvents4 = max(finalCombinedDataTempPA(:,4));
    

    if maxNumEvents4 > 0
    TempCellArrayEventPA4=[];      
    TempCellArrayEventPA4{maxNumEvents4,13}=[];     
        for i = 1: maxNumEvents4
        locs3 = find(finalCombinedDataTempPA(:,4) == i);
         if numel(locs3) > 5
            FiberDataArrayPA.NonActive{i} = partialAnalysisData(locs3, 4);
            finalCombinedDataTempPA4 = partialAnalysisData(locs3, 6);
            [pksPAEvent3,locsPAEvent4]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);%May23,2023
            TempCellArrayEventPA4{i,1}=basename1;
            TempCellArrayEventPA4{i,2}=basename2;
            TempCellArrayEventPA4{i,3}="OtherZone";
            TempCellArrayEventPA4{i,4}=mean(partialAnalysisData(locs3, 6));%May23,2023
            TempCellArrayEventPA4{i,5}=mean(sigmaplotF(locs3));% from line 802  % mean Z-score 
            TempCellArrayEventPA4{i,6}=numel(locsPAEvent4);
            TempCellArrayEventPA4{i,7}=mean(finalCombinedDataTempPA4(locsPAEvent4));
            TempCellArrayEventPA4{i,8} =str2num(reply1);
            TempCellArrayEventPA4{i,9}=numel(finalCombinedDataTempPA4);
            
            sigmaplotFtemp4=sigmaplotF(locs3); %HT 02/02/2024
            TempCellArrayEventPA4{i,10}=mean(sigmaplotFtemp4(locsPAEvent4)); % Mean Peak Height Z-Score

            [pksPAEvent04,locsPAEvent04]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %02/27/2025 PKN
            TempCellArrayEventPA4{i,11}= numel(locsPAEvent04);%02/27/2025 PKN
            TempCellArrayEventPA4{i,12}= mean(finalCombinedDataTempPA4(locsPAEvent04));%03/03/2025 PKN
            TempCellArrayEventPA4{i,13}=mean(sigmaplotFtemp4(locsPAEvent04));%03/03/2025 PKN

         end
        end
    else 
    TempCellArrayEventPA4=[];     
    end
    
    
    
    
    
    
    
    
   TempCellArrayAllEventsPA=[];
   TempCellArrayAllEventsPA = [TempCellArrayEventPA1 ; TempCellArrayEventPA2; TempCellArrayEventPA3 ; TempCellArrayEventPA4];

   TempCellLabelEventsPA={"basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z", "numSig" "meanSigDFF" "MeanSigDFF in Z"};
   TempCellLabelEventsPA1=["basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z", "numSig" "meanSigDFF" "MeanSigDFF in Z"];


   TableTempCellLabelEventsPA=cell2table(TempCellLabelEventsPA, "VariableNames", TempCellLabelEventsPA1);
   datafilenameEventsPA=['125_d-SLR_SAMPLE_events_partial_081325' '_' num2str(startpoint) '_' num2str(endpoint) '.csv'];
    if isfile(datafilenameEventsPA);
      
    else
    writetable(TableTempCellLabelEventsPA,datafilenameEventsPA);
    end
    
     CellTableEventsPA=table2cell(readtable(datafilenameEventsPA));
     CellTableNewEventsPA=[CellTableEventsPA;TempCellArrayAllEventsPA]; %%%NOTE! May23,2023
     TableNewEventsPA=cell2table(CellTableNewEventsPA, "VariableNames", TempCellLabelEventsPA1);
     writetable(TableNewEventsPA,datafilenameEventsPA);


    datafilenameEventsPAMat = [basename1 'PA_events.mat'];

    %matlab file for all 3 events in a data structure
    save(datafilenameEventsPAMat, "FiberDataArrayPA");
    
    
   % calculate average fluorescence during the zone 
    
    TotalZone1CountF=sum(partialAnalysisData(:,7));
    TempArrayF=zeros(TotalZone1CountF,1);
    
    TotalZone2CountF=sum(partialAnalysisData(:,8));
    TempArray2F=zeros(TotalZone2CountF,1);
    
    TotalZone3CountF=sum(partialAnalysisData(:,9));
    TempArray3F=zeros(TotalZone3CountF,1);
        
    TotalOtherZoneCountF=sum(partialAnalysisData(:,10));
    TempArray4F=zeros(TotalOtherZoneCountF,1);
    
    
    
    nn=1; %counter 
    nnn=1; %counter 
    nnnn=1; %counter 
    nnnnn=1; %11/17/2023
    for iii=1:lengthF 
     if partialAnalysisData(iii,7)==1;
         TempArrayF(nn,1)=partialAnalysisData(iii,6);
         nn=nn+1;
     elseif partialAnalysisData(iii,8)==1
         TempArray2F(nnn,1)=partialAnalysisData(iii,6);
         nnn=nnn+1;
     elseif partialAnalysisData(iii,9)==1
         TempArray3F(nnnn,1)=partialAnalysisData(iii,6);
         nnnn=nnnn+1;
     else
         TempArray4F(nnnnn,1)=partialAnalysisData(iii,6);
         nnnnn=nnnnn+1;
     end
    end
    


    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    nnnnn=1; %11/17
    for iii=1:numpks2;
        NNN=peaktable2(iii,1);
        if partialAnalysisData(NNN,7)==1;
            Zone1PeakF(nn,1)=peaktable2(iii,3);
            nn=nn+1;
        elseif partialAnalysisData(NNN,8)==1;
            Zone2PeakF(nnn,1)=peaktable2(iii,3);
            nnn=nnn+1;
        elseif partialAnalysisData(NNN,9)==1
            Zone3PeakF(nnnn,1)=peaktable2(iii,3);
            nnnn=nnnn+1;
        else
            OtherZonePeakF(nnnnn,1)=peaktable2(iii,3);
            nnnnn=nnnnn+1;
        end 
    end
   
   if nn==1;
       Zone1PeakF=[];
   end
   
   if nnn==1;
       Zone2PeakF=[];
   else
   end
   
   if nnnn==1;
       Zone3PeakF=[];
   else
   end
   
   if nnnnn==1;
       OtherZonePeakF=[];
   else
   end
    
   DFFZone1Mean=mean(TempArrayF);
   DFFZone1Stdev=std(TempArrayF);
   DFFZone1Counts=size(TempArrayF,1);
 
   DFFZone2Mean=mean(TempArray2F);
   DFFZone2Stdev=std(TempArray2F);
   DFFZone2Counts=size(TempArray2F,1); 
    
        
   DFFZone3Mean=mean(TempArray3F);
   DFFZone3Stdev=std(TempArray3F);
   DFFZone3Counts=size(TempArray3F,1); 
   
   DFFOtherZoneMean=mean(TempArray4F);
   DFFOtherZoneStdev=std(TempArray4F);
   DFFOtherZoneCounts=size(TempArray4F,1); 
   
   
   
    
   Zone1PeakCounts=size(Zone1PeakF,1);
   NormalizedZone1PeakCounts=Zone1PeakCounts/DFFZone1Counts*20;
   Zone1PeakHeight=mean(Zone1PeakF);
    
   Zone2PeakCounts=size(Zone2PeakF,1);
   NormalizedZone2PeakCounts=Zone2PeakCounts/DFFZone2Counts*20;
   Zone2PeakHeight=mean(Zone2PeakF);
    
   Zone3PeakCounts=size(Zone3PeakF,1);
   NormalizedZone3PeakCounts=Zone3PeakCounts/DFFZone3Counts*20;
   Zone3PeakHeight=mean(Zone3PeakF);
   
   
   OtherZonePeakCounts=size(OtherZonePeakF,1);
   NormalizedOtherZonePeakCounts=OtherZonePeakCounts/DFFOtherZoneCounts*20;
   OtherZonePeakHeight=mean(OtherZonePeakF);  
   
   
   
    
   TempCellArray{1,30}=[]; 
   TempCellArray{1,1}=basename1;
   TempCellArray{1,2}=basename2;
   
   TempCellArray{1,3}=NSegment;
   
   TempCellArray{1,4}=startpoint;
   TempCellArray{1,5}=endpoint;
  
   TempCellArray{1,6}=DFFZone1Mean;
   TempCellArray{1,7}=DFFZone1Stdev;
   TempCellArray{1,8}=DFFZone1Counts;
   
   TempCellArray{1,9}=DFFZone2Mean;
   TempCellArray{1,10}=DFFZone2Stdev;
   TempCellArray{1,11}=DFFZone2Counts;
   
   TempCellArray{1,12}=DFFZone3Mean;
   TempCellArray{1,13}=DFFZone3Stdev;
   TempCellArray{1,14}=DFFZone3Counts;
   

   TempCellArray{1,15}=DFFOtherZoneMean;
   TempCellArray{1,16}=DFFOtherZoneStdev;
   TempCellArray{1,17}=DFFOtherZoneCounts;
   
   
   TempCellArray{1,18}=Zone1PeakCounts;
   TempCellArray{1,19}=NormalizedZone1PeakCounts;
   TempCellArray{1,20}=Zone1PeakHeight;
   
   TempCellArray{1,21}=Zone2PeakCounts;
   TempCellArray{1,22}=NormalizedZone2PeakCounts;
   TempCellArray{1,23}=Zone2PeakHeight;
   
   TempCellArray{1,24}=Zone3PeakCounts;
   TempCellArray{1,25}=NormalizedZone3PeakCounts;
   TempCellArray{1,26}=Zone3PeakHeight;
   
   TempCellArray{1,27}=OtherZonePeakCounts;
   TempCellArray{1,28}=NormalizedOtherZonePeakCounts;
   TempCellArray{1,29}=OtherZonePeakHeight;
   
   
   
   TempCellArray{1,30}=reply1;

  
    
   % TempCellLabel={'basename1' 'basename' 'NSegment' 'datastart' 'dataend' 'DFFZoneMean' 'DFFZoneStdev' 'DFFZoneCounts' 'DFFCornerZoneMean' 'DFFCornerZoneStdev' 'DFFCornerZoneCounts' 'DFFNonZoneMean' 'DFFNonZoneStdev' 'DFFNonZoneCounts' 'ZonePeakCounts' 'NormalizedZonePeakCounts' 'ZonePeakHeight' 'CornerZonePeakCounts' 'NormalizedCornerZonePeakCounts' 'CornerZonePeakHeight' 'NonZonePeakCounts' 'NormalizedNonZonePeakCounts' 'NonZonePeakHeight'};
%    TempCellLabel={"basename1" "basename" "NSegment" "datastart" "dataend" "DFFZoneMean" "DFFZoneStdev" "DFFZoneCounts" "DFFCornerZoneMean" "DFFCornerZoneStdev" "DFFCornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneStdev" "DFFNonZoneCounts" "ZonePeakCounts" "NormalizedZonePeakCounts" "ZonePeakHeight" "CornerZonePeakCounts" "NormalizedCornerZonePeakCounts" "CornerZonePeakHeight" "NonZonePeakCounts" "NormalizedNonZonePeakCounts" "NonZonePeakHeight" "TypeOfAnalysisControl"};
%    TempCellLabel1=["basename1" "basename" "NSegment" "datastart" "dataend" "DFFZoneMean" "DFFZoneStdev" "DFFZoneCounts" "DFFCornerZoneMean" "DFFCornerZoneStdev" "DFFCornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneStdev" "DFFNonZoneCounts" "ZonePeakCounts" "NormalizedZonePeakCounts" "ZonePeakHeight" "CornerZonePeakCounts" "NormalizedCornerZonePeakCounts" "CornerZonePeakHeight" "NonZonePeakCounts" "NormalizedNonZonePeakCounts" "NonZonePeakHeight" "TypeOfAnalysisControl"];

 TempCellLabel={"basename1" "basename" "NSegment" "datastart" "dataend" "DFFZone1Mean" "DFFZone1Stdev" "DFFZone1Counts" "DFFZone2Mean" "DFFZone2Stdev" "DFFZone2Counts" "DFFZone3Mean" "DFFZone3Stdev" "DFFZone3Counts" "DFFOtherZoneMean" "DFFOtherZoneStdev" "DFFOtherZoneCounts" "Zone1PeakCounts" "NormalizedZone1PeakCounts" "Zone1PeakHeight" "Zone2PeakCounts" "NormalizedZone2PeakCounts" "Zone2PeakHeight" "Zone3PeakCounts" "NormalizedZone3PeakCounts" "Zone3PeakHeight"  "OtherZonePeakCounts" "NormalizedOtherZonePeakCounts" "OtherZonePeakHeight" "TypeOfAnalysisControl"};


TempCellLabel1=["basename1" "basename" "NSegment" "datastart" "dataend" "DFFZone1Mean" "DFFZone1Stdev" "DFFZone1Counts" "DFFZone2Mean" "DFFZone2Stdev" "DFFZone2Counts" "DFFZone3Mean" "DFFZone3Stdev" "DFFZone3Counts" "DFFOtherZoneMean" "DFFOtherZoneStdev" "DFFOtherZoneCounts" "Zone1PeakCounts" "NormalizedZone1PeakCounts" "Zone1PeakHeight" "Zone2PeakCounts" "NormalizedZone2PeakCounts" "Zone2PeakHeight" "Zone3PeakCounts" "NormalizedZone3PeakCounts" "Zone3PeakHeight"  "OtherZonePeakCounts" "NormalizedOtherZonePeakCounts" "OtherZonePeakHeight" "TypeOfAnalysisControl"];
  
 
  %    datafilename=['FiberDataExport4' '.csv'];
%     if ~isempty(datafilename);
%         xlswrite(datafilename,TempCellLabel,'Sheet1');
%     else
%     end
% 
%     [success,message]=xlsappend(datafilename,TempCellArray);
   TableTempCellLabel=cell2table(TempCellLabel, "VariableNames", TempCellLabel1);
   datafilename=['125_d-SLR_SAMPLE_summary_PARTIAL_081325' '_' num2str(startpoint) '_' num2str(endpoint) '.csv'];
    if isfile(datafilename);
      
    else
    writetable(TableTempCellLabel,datafilename);
    end
    
     CellTable=table2cell(readtable(datafilename));
 
   
     writetable(TableNew,datafilename);
 else
 end   
     
   
   
 
    
 %%   %
function see = myfitFPSY(params,Iso,Calcium)

%%%%% Iso is 410nm Signal 
%%%%% Calcium  is 460nm Signal 

a = params(1);
b = params(2);

Fitted_Curve = a*Iso+b;
Error_Vector = Fitted_Curve - Calcium;
see = sum(Error_Vector.^2);

end
% end

function[Y]=CountConsecutiveOnes(X);
%%%%%% [ 0 0 0 1 1 0 0 0 1 1 1 1 0 0 1 1 1 ] 
%%%%%% [ 0 0 0 2 2 0 0 0 4 4 4 4 0 0 2 2 1 ]
%%%%%% Note: Last number remain unchanged becasue it will calculate one
%%%%%% position ahead. Just like DIFF signal, the new array is one number less than
%%%%%% original - but keep the data length unchanged for convenience. 
%%%%%%  Feb 24, 2017, HAJI 



lnX=length(X);
Y=X; %%%% New Array 
C=zeros(lnX); %%%% Investigation Check  

for i=1:lnX
    if (X(i)==1 && C(i)==0);
        
        m=i+1;
        count=1;
        if m < lnX;
        while (X(m)==1 && m < lnX)
            count=count+1;
            m=m+1;
        end
        for k=i:m-1;
        Y(k)=count;
        C(k)=1;
        end
        else
            Y(i)=1;
            C(i)=1;
        end
    else
        C(i)=1;
    end
end

end