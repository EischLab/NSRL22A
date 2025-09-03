   
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
%%

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

%select a video data
[filename2,pathname2]=uigetfile('*.csv');
basename2=strrep(filename2,'.csv','');
addpath pathname2;
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

NSegment=1;
startpoint=1;
endpoint=lengthFiberData;
    
    vlength=size(videoDataRaw,1);
    videoDataClean=zeros(vlength,8);%% E1=1,2,3, E2=0.033,0.066 in sec, E3=original zone binary , E4 reversed binary E5,processed 
    videoDataClean(:,1)=[1:vlength];
    videoDataClean(:,2)=0.03333*(videoDataClean(:,1)-1);
    videoDataClean(:,4)=1;
    
    
% identify interaction zone     
    for iii=1:vlength;
        if videoDataRaw{iii,7}=="True";
             videoDataClean(iii,3)=1; % zone
             videoDataClean(iii,4)=0; % reverse (non zone) for cleaning 
        elseif videoDataRaw{iii,7}=="TRUE";
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
    
    
   %corner zone 
   for iii=1:vlength;
        if videoDataRaw{iii,5}=="True" | videoDataRaw{iii,6}=="True";
             videoDataClean(iii,7)=1; % corner zones
         elseif videoDataRaw{iii,5}=="TRUE" | videoDataRaw{iii,6}=="TRUE";
             videoDataClean(iii,7)=1; % corner zones
        end
   end
   %other area
   for iii=1:vlength;
        if videoDataClean(iii,6)==0 && videoDataClean(iii,7)==0;
             videoDataClean(iii,8)=1; % other area
         else
        end
    end
    
 %   
 
   

 %Run for different control for peak detection +  
 % Run to create full figure analysis

 reply1 = input(['Do you want normal peak detection using original 415nm data (1),' ...
     'using smoothed 415 nm data (2), or using only 470nm data (3):'],"s");
  
    
 % Fiberphotometry data (fluorescecne@470nm)  - baseline correction by fluorescecne@410nm.  
 % initilize array
 
    Fitted=[];
    DeltaFoverF=[];
    combinedFibVidData=zeros(lengthFiberData,10);% output array 
    fiberZScores=zeros(lengthFiberData,2);
    finalCombinedData=zeros(lengthFiberData,7);
    
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
    xlabel('Frame Number');
    ylabel('Wavelength');
    plot(tempControl,'Color','blue');
    hold on
    plot(Raw,'Color','green');
    hold on
    plot(Control ,'Color','magenta');
    
    
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

    %%
  
    
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
    xlabel('Frame Number');
    ylabel('470nm');
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
        
        combinedFibVidData(jj,8)=videoDataClean(vindex,6);
        combinedFibVidData(jj,9)=videoDataClean(vindex,7);
        combinedFibVidData(jj,10)=videoDataClean(vindex,8);
    end
    
    
    finalCombinedData(:,5:7)=combinedFibVidData(:,8:10);  % interaction zone, corner zone , other 
    
    
  
     
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
    [pks0,locs0]=findpeaks(DeltaFoverF-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);% 02/28/25 PKN
    
    %[pks1,locs1]=findpeaks(sigmaplot,'MinPeakHeight',2,'MinPeakDistance',3);
    numpks=length(locs1);
    numpks0=length(locs0);%02/27/2025 PKN
    
    
    peaktable=zeros(numpks,3);%initialize peaktable 
    peaktable(:,1)=locs1;
    peaktable(:,2)=pks1;
    
    dff0pks1=DeltaFoverF(locs1);
    peaktable(:,3)=dff0pks1;
    peaktable(:,4) = fiberZScores(locs1,2);
    peaktable(:,5) = combinedFibVidData(locs1,3);
    locs2 = combinedFibVidData(locs1,3); % extracts corresponding time for an occurring peak

    %initialize a table containing all signals before threshold %02/27/2025 PKN
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
    xlabel('Frame Number');
    ylabel('Z-Score(DFF)');
    plot(fiberZScores(:,1),fiberZScores(:,2),'Color','black');
    hold on 
    for ii=1:numpks;
        scatter(locs1(ii),pks1(ii),5,'blue');
        hold on
    end
    
    
    figure(H0);
    subplot(5,2,[5,6])
    plot(combinedFibVidData(:,1),combinedFibVidData(:,6),'Color','black');
    xlabel('Frame Number');
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
    
 
    %second graph with time in s as x-axis (time from FP) %12/2024 PKN

    figure(H0);
    subplot(5,2,[7,8])
    title('Z-score after threshold')
    xlabel('Time in seconds');

   plot(combinedFibVidData(:,3),combinedFibVidData(:,6),'Color','black');
    hold on 
    for ii=1:numpks;
        scatter(locs2(ii),dff0pks1(ii),5,'blue');
        hold on
    end
    %08/05/25_PKN


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
   
    

% ploting a graph showing all signals before threshold with the zone
% 02/28/25PKN
    figure(H0);
    subplot(5,2,[9,10])
    title('Z-score before threshold')
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
    xlim([0, max(combinedFibVidData(:,1))]);
  
   figurefilename2=[basename1 '-.png'];
   saveas(H0,figurefilename2);

   
   %creating separate events for averaging

    finalCombinedDataTemp(:,1) = bwlabel(logical(finalCombinedData(:,5)));
    maxNumEvents1 = max(finalCombinedDataTemp(:,1));
    
   

    if maxNumEvents1 ~=0
        
    TempCellArrayEvent1{maxNumEvents1,10}=[]; 
    
        for i = 1: maxNumEvents1
        locs3 = find(finalCombinedDataTemp(:,1) == i);
        if numel(locs3) > 5;
        FiberDataArray.Active{i} = finalCombinedData(locs3, 4);
        ZFiberDataArray.Active{i} = fiberZScores(locs3, 2);
        
        finalCombinedDataTemp1 = finalCombinedData(locs3, 4);
        [pksEvent1,locsEvent1]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
        TempCellArrayEvent1{i,1}=basename1;
        TempCellArrayEvent1{i,2}=basename2;
        TempCellArrayEvent1{i,3}="Active";
        TempCellArrayEvent1{i,4}=mean(finalCombinedData(locs3, 4));
        TempCellArrayEvent1{i,5}=mean(sigmaplot(locs3));
        TempCellArrayEvent1{i,6}=numel(locsEvent1);
        TempCellArrayEvent1{i,7}=mean(finalCombinedDataTemp1(locsEvent1));
        TempCellArrayEvent1{i,8} = str2num(reply1);
        TempCellArrayEvent1{i,9} =numel(finalCombinedDataTemp1);
        sigmatemp1=sigmaplot(locs3); %02/02/2024 HT
        TempCellArrayEvent1{i,10}=mean(sigmatemp1(locsEvent1));%02/02/2024 HT
        TempCellArrayEvent1{i,11}= numel(sigmatemp1(locsEvent1));
    
       [pksEvent01,locsEvent01]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %02/27/2025 PKN
       TempCellArrayEvent1{i,12}= numel(locsEvent01);%02/27/2025 PKN
       TempCellArrayEvent1{i,13}= mean(finalCombinedDataTemp1(locsEvent01));%03/03/2025 PKN
       TempCellArrayEvent1{i,14}=mean(sigmatemp1(locsEvent01));%03/03/2025 PKN
        end
        end
    else 
        TempCellArrayEvent1=[];
    end


    finalCombinedDataTemp(:,2) = bwlabel(logical(finalCombinedData(:,6)));
    maxNumEvents2 = max(finalCombinedDataTemp(:,2));
   

    if maxNumEvents2 ~=0
        
     TempCellArrayEvent2{maxNumEvents2,10}=[]; 
     
        for i = 1: maxNumEvents2
        locs3 = find(finalCombinedDataTemp(:,2) == i);
        if numel(locs3) > 5
        FiberDataArray.Corner{i} = finalCombinedData(locs3, 4);
        ZFiberDataArray.Corner{i} = fiberZScores(locs3, 2);
        finalCombinedDataTemp2 = finalCombinedData(locs3, 4);
        
        [pksPAEvent2,locsEvent2]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
        TempCellArrayEvent2{i,1}=basename1;
        TempCellArrayEvent2{i,2}=basename2;
        TempCellArrayEvent2{i,3}="Corner";
        TempCellArrayEvent2{i,4}=mean(finalCombinedData(locs3, 4));
        TempCellArrayEvent2{i,5}=mean(sigmaplot(locs3));
        TempCellArrayEvent2{i,6}=numel(locsEvent2);
        TempCellArrayEvent2{i,7}=mean(finalCombinedDataTemp2(locsEvent2));
        TempCellArrayEvent2{i,8} = str2num(reply1);
        TempCellArrayEvent2{i,9} = numel(finalCombinedDataTemp2);
      
        sigmatemp2=sigmaplot(locs3); %02/02/2024 HT
        TempCellArrayEvent2{i,10}=mean(sigmatemp2(locsEvent2)); %02/02/2024 HT
        TempCellArrayEvent2{i,11}= numel(sigmatemp2(locsEvent2));
        [pksEvent02,locsEvent02]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%02/27/2025 PKN
        TempCellArrayEvent2{i,12}= numel(finalCombinedData(locsEvent02, 4));%02/27/2025 PKN
        TempCellArrayEvent2{i,13}= mean(finalCombinedData(locsEvent02, 4));%03/03/2025 PKN
        TempCellArrayEvent2{i,14}=mean(sigmatemp2(locsEvent02));%03/03/2025 PKN
        end
        end
    else 
    TempCellArrayEvent2=[];     
    end

    finalCombinedDataTemp(:,3) = bwlabel(logical(finalCombinedData(:,7)));
    maxNumEvents3 = max(finalCombinedDataTemp(:,3));
    

    if maxNumEvents3 ~=0
        
    TempCellArrayEvent3{maxNumEvents3,10}=[]; 
    
        for i = 1: maxNumEvents3
        locs3 = find(finalCombinedDataTemp(:,3) == i);
        FiberDataArray.NonActive{i} = finalCombinedData(locs3, 4);
        ZFiberDataArray.NonActive{i} = fiberZScores(locs3, 2);
        if numel(locs3) > 5
            finalCombinedDataTemp3 = finalCombinedData(locs3, 4);
            [pksPAEvent3,locsEvent3]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
            TempCellArrayEvent3{i,1}=basename1;
            TempCellArrayEvent3{i,2}=basename2;
            TempCellArrayEvent3{i,3}="NonActive";
            TempCellArrayEvent3{i,4}=mean(finalCombinedData(locs3, 4));
            TempCellArrayEvent3{i,5}=mean(sigmaplot(locs3));
            TempCellArrayEvent3{i,6}=numel(locsEvent3);
            TempCellArrayEvent3{i,7}=mean(finalCombinedDataTemp3(locsEvent3));
            TempCellArrayEvent3{i,8} = str2num(reply1);
            TempCellArrayEvent3{i,9} =numel(finalCombinedDataTemp3);
            
            sigmatemp3=sigmaplot(locs3); %02/02/2024 HT
            TempCellArrayEvent3{i,10}=mean(sigmatemp3(locsEvent3));
            TempCellArrayEvent3{i,11}= numel(sigmatemp3(locsEvent3));
            [pksEvent03,locsEvent03]=findpeaks(finalCombinedData(locs3, 4)-meanF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%02/27/2025 PKN
            TempCellArrayEvent3{i,12}= numel(finalCombinedData(locsEvent03, 4));%02/27/2025 PKN
            TempCellArrayEvent3{i,13}= mean(finalCombinedData(locsEvent03, 4));%03/03/2025 PKN
            TempCellArrayEvent3{i,14}=mean(sigmatemp3(locsEvent03));%03/03/2025 PKN

        end
        end
    else 
        
    TempCellArrayEvent3=[];    
    end

    TempCellArrayAllEvents = [TempCellArrayEvent1 ; TempCellArrayEvent2; TempCellArrayEvent3];

   TempCellLabelEvents={"basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z", "numofSig(zs)" "numSig" "meanSigDFF" "MeanSigDFF in Z"};
   TempCellLabelEvents1=["basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z" "numofSig(zs)" "numSig" "meanSigDFF" "MeanSigDFF in Z"];


   TableTempCellLabelEvents=cell2table(TempCellLabelEvents, "VariableNames", TempCellLabelEvents1);
   datafilenameEvents=['125_xs-SLR_TEST_events_all_080525' '.csv'];
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
%%
%{

    %--------- Plotting Fluorescence Activity across events while in Active Zone---------------------------------------
    H1=figure('Name',basename2);
    figure(H1);
    
    subplot(2,1,1);
    hold on; % Allows multiple plots on the same figure

   for i = 1:length(ZFiberDataArray.Active)
    plot(ZFiberDataArray.Active{i}, 'LineWidth', 1.5);
  end
   
   ylabel('z-score');
   title('Activity Across Events in Interaction Zone');
   
   grid on;
   hold off;

   subplot(2,1,2);
   hold on;
   for i = 1:length(FiberDataArray.Active)
    plot(FiberDataArray.Active{i}, 'LineWidth', 1.5);
  end
   
   ylabel('Df/F');
   
   grid on;
   hold off;

   figurefilename1=[basename1,'InteractionZone_ALL_DF_F.png'];
   saveas(H1,figurefilename1);
  
   
   % ----------Plotting Fluorescence Activity across events while Corner Zone-------------------------------------- -
    H2=figure('Name',basename2);
    figure(H2);
    
    subplot(2,1,1);
    hold on; % Allows multiple plots on the same figure
   for i = 1:length(ZFiberDataArray.Corner)
    plot(ZFiberDataArray.Corner{i}, 'LineWidth', 1.5);
  end
   
   ylabel('z-score');
   title('Activity Across Events in Corner Zone');
   
   grid on;
   hold off;
   
   subplot(2,1,2);
   hold on;
   for i = 1:length(FiberDataArray.Corner)
    plot(FiberDataArray.Corner{i}, 'LineWidth', 1.5);
  end
   
   ylabel('Df/F');
   
   grid on;
   hold off;
   
   
   figurefilename3=[basename1,'CornerZone_ALL_DF_f.png'];
   saveas(H2,figurefilename3);

%-----------Creating one figure with mulitple events- Interaction Zone %030625 PKN-----------------------------------

    H3=figure('Name',basename2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    figure(H3);
    numEvents = length(ZFiberDataArray.Active); % Number of events
    numCols = 5; % Set the number of columns (change to 3 if needed)
    numRows = ceil(numEvents / numCols); % Compute required rows

    for i = 1:length(ZFiberDataArray.Active)
    ymin = min(ZFiberDataArray.Active{i});
end

    xMin = 0;
    xMax = max(cellfun(@length, ZFiberDataArray.Active)); % Longest event
    ymin = min(ZFiberDataArray.Active{i});
    ymax = max(ZFiberDataArray.Active{i});
    ylim([ymin ymax]);
    
validEventIdx = 1; % Counter for valid plots
for i = 1:numEvents
    data = ZFiberDataArray.Active{i}; % Extract event data

    % Skip if data is empty or all values are NaN
    if isempty(data) || all(isnan(data))
        continue;
    end
    subplot(numRows, numCols, i);
    plot(ZFiberDataArray.Active{i}, 'b-', 'LineWidth', 1);
   
    title(['Event Number ', num2str(i)]);
    grid on; 
end

sgtitle('Fluorescence Activity Across Events In Interaction Zone (Zscore)'); % Super title
figurefilename4=[basename1,'InteractionZone_Individual_Zscore.png'];
   saveas(H3,figurefilename4);

%Creating one figure with mulitple events- Corner Zone %031125 PKN
    H4=figure('Name',basename2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    figure(H4);
  
    numEvents = length(ZFiberDataArray.Corner); % Number of events
    numCols = 5; % Set the number of columns (change to 3 if needed)
    numRows = ceil(numEvents / numCols); % Compute required rows
    validEventIdx = 1; % Counter for valid plots
for i = 1:numEvents
    data = ZFiberDataArray.Corner{i}; % Extract event data

    % Skip if data is empty or all values are NaN
    if isempty(data) || all(isnan(data))
        continue;
    end
    xMin = 0;
    xMax = max(cellfun(@length, ZFiberDataArray.Corner)); % Longest event
    ymin = min(ZFiberDataArray.Corner{i});
    ymax = max(ZFiberDataArray.Corner{i});
    ylim([ymin ymax]);
    

    subplot(numRows, numCols, i);
    plot(ZFiberDataArray.Corner{i}, 'b-', 'LineWidth', 1);
   
    title(['Event Number ', num2str(i)]);
    grid on;
end

sgtitle('Fluorescence Activity Across Events In Corner Zone (Zscore)'); % Super title
figurefilename5=[basename1,'CornerZone_Individual_Zscore.png'];
   saveas(H4,figurefilename5);
   % This is the "next section" of your script.
% Code here will execute regardless of the user's choice.
% For example:
disp('Continuing with the rest of the script...');
%}
   %%
%calculate average fluorescence during the zone 
    
    TotalActiveCount=sum(finalCombinedData(:,5));
    TempArray=zeros(TotalActiveCount,1);
    
    TotalCornerCount=sum(finalCombinedData(:,6));
    TempArray2=zeros(TotalCornerCount,1);
    
    TotalNonActiveCount=sum(finalCombinedData(:,7));
    TempArray3=zeros(TotalNonActiveCount,1);

    
    nn=1; %counter 
    nnn=1; %counter 
    nnnn=1; %counter 
    for iii=1:lengthFiberData 
     if finalCombinedData(iii,5)==1;
         TempArray(nn,1)=finalCombinedData(iii,4);
         nn=nn+1;
     elseif finalCombinedData(iii,6)==1
         TempArray2(nnn,1)=finalCombinedData(iii,4);
         nnn=nnn+1;
     else
         TempArray3(nnnn,1)=finalCombinedData(iii,4);
         nnnn=nnnn+1;
     end
    end

    
    %Finding DF/F signal amplitude (before thresold) %03/03/25 PKN 
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks0;
        NNN=Sigtable(iii,1);
        if finalCombinedData(NNN,5)==1;
            ActiveZoneSig(nn,1)=Sigtable(iii,3);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            CornerZoneSig(nnn,1)=Sigtable(iii,3);
            nnn=nnn+1;
        else
            NonActiveZoneSig(nnnn,1)=Sigtable(iii,3);
            nnnn=nnnn+1;
        end 
    end
        
  %Finding Z-score of DF/F signal amplitude (before thresold) %03/03/25 PKN 
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks0;
        NNN=Sigtable(iii,1);
        if finalCombinedData(NNN,5)==1;
            ZscoreActiveZoneSig(nn,1)=Sigtable(iii,4);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            ZscoreCornerZoneSig(nnn,1)=Sigtable(iii,4);
            nnn=nnn+1;
        else
            ZscoreNonActiveZoneSig(nnnn,1)=Sigtable(iii,4);
            nnnn=nnnn+1;
        end 
    end
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks;
        NNN=peaktable(iii,1);
        if finalCombinedData(NNN,5)==1;
            ActiveZonePeak(nn,1)=peaktable(iii,3);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            CornerZonePeak(nnn,1)=peaktable(iii,3);
            nnn=nnn+1;
        else
            NonActiveZonePeak(nnnn,1)=peaktable(iii,3);
            nnnn=nnnn+1;
        end 
    end
% Finding Z-score of peak amplitude  %03/03/25 PKN     
  nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;

    for iii=1:numpks;
        NNN=peaktable(iii,1);
        if finalCombinedData(NNN,5)==1;
            ZscoreActiveZonePeak(nn,1)=peaktable(iii,4);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            ZscoreCornerZonePeak(nnn,1)=peaktable(iii,4);
            nnn=nnn+1;
        else
            ZscoreNonActiveZonePeak(nnnn,1)=peaktable(iii,4);
            nnnn=nnnn+1;
        end 
    end  
   if nn==1;
       ActiveZonePeak=[];
   end
   
   if nnn==1;
       CornerZonePeak=[];
   else
   end
   
   if nnnn==1
       NonActiveZonePeak=[];
   else
   end
   if nn==1;
       ZscoreActiveZonePeak=[];
   end
   
   if nnn==1;
       ZscoreCornerZonePeak=[];
   else
   end
   
   if nnnn==1
       ZscoreNonActiveZonePeak=[];
   else
   end
    
   ActiveZoneSig = []; 

   if ~isempty(ActiveZoneSig)
    % The array has data, so calculate the mean.
    DFFZoneSigMean = mean(ActiveZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
    DFFZoneSigMean = NaN; 
   end

   ZscoreActiveZoneSig = []; 

   if ~isempty(ActiveZoneSig)
    % The array has data, so calculate the mean.
    ZscoreDFFZoneSigMean = mean(ZscoreActiveZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
    ZscoreDFFZoneSigMean = NaN; 
    end

    CornerZoneSig=[]
    if ~isempty(CornerZoneSig)
    % The array has data, so calculate the mean.
      DFFCornerZoneSigMean= mean(CornerZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
      DFFCornerZoneSigMean = NaN; 
    end

    ZscoreCornerZoneSig=[]
    if ~isempty(ZscoreCornerZoneSig)
    % The array has data, so calculate the mean.
    ZscoreDFFCornerZoneSigMean= mean(ZscoreCornerZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
    ZscoreDFFCornerZoneSigMean = NaN; 
    end

   DFFZoneMean=mean(TempArray);
   DFFZoneSigMean=mean(ActiveZoneSig);
   ZscoreDFFZoneSigMean=mean(ZscoreActiveZoneSig);
   DFFZoneStdev=std(TempArray);
   ZoneCounts=size(TempArray,1);
   DFFZoneSigCounts=size(ActiveZoneSig,1);
   ZoneSignalRate=DFFZoneSigCounts/ZoneCounts*20;% peak#/sec
  
   DFFCornerZoneMean=mean(TempArray2);
   DFFCornerZoneSigMean=mean(CornerZoneSig);
   ZscoreDFFCornerZoneSigMean= mean(ZscoreCornerZoneSig);
   DFFCornerZoneStdev=std(TempArray2);
   CornerZoneCounts=size(TempArray2,1);
   DFFCornerZoneSigCounts=size(CornerZoneSig,1);
   CornerZoneSignalRate=DFFCornerZoneSigCounts/CornerZoneCounts*20;% peak#/sec  
        
   DFFNonZoneMean=mean(TempArray3);
   DFFNonZoneSigMean=mean(NonActiveZoneSig);
   ZscoreDFFNonZoneSigMean=mean(ZscoreNonActiveZoneSig);
   DFFNonZoneStdev=std(TempArray3);
   NonZoneCounts=size(TempArray3,1); 
   DFFNonZoneSigCounts=size(NonActiveZoneSig,1);
   NonZoneSignalRate=DFFNonZoneSigCounts/NonZoneCounts*20;% peak#/sec
  
    
   ZonePeakCounts=size(ActiveZonePeak,1);
   RateZonePeakCounts=ZonePeakCounts/ZoneCounts*20;% peak#/sec
   ZonePeakHeight=mean(ActiveZonePeak);
   ZscoreZonePeakHeight=mean(ZscoreActiveZonePeak);
    
   CornerPeakCounts=size(CornerZonePeak,1);
   RateCornerPeakCounts=CornerPeakCounts/CornerZoneCounts*20;
   CornerPeakHeight=mean(CornerZonePeak);

   ZscoreCornerZonePeakHeight=mean(ZscoreCornerZonePeak);
    
   NonZonePeakCounts=size(NonActiveZonePeak,1);
   RateNonZonePeakCounts=NonZonePeakCounts/NonZoneCounts*20;
   NonZonePeakHeight=mean(NonActiveZonePeak);
   ZscoreNonZonePeakHeight=mean(ZscoreNonActiveZonePeak); 
   
    
   TempCellArray{1,36}=[]; 
   TempCellArray{1,1}=basename1;
   TempCellArray{1,2}=basename2;
   
   TempCellArray{1,3}=NSegment;
   
   TempCellArray{1,4}=startpoint;
   TempCellArray{1,5}=endpoint;
  
   TempCellArray{1,6}=DFFZoneMean;
   TempCellArray{1,7}=DFFZoneSigMean;
   TempCellArray{1,8}=ZscoreDFFZoneSigMean;
   TempCellArray{1,9}=DFFZoneSigCounts;
   TempCellArray{1,10}=ZoneSignalRate;
   TempCellArray{1,11}=DFFZoneStdev;
   TempCellArray{1,12}=ZoneCounts;
   
   TempCellArray{1,13}=DFFCornerZoneMean;
   TempCellArray{1,14}=DFFCornerZoneSigMean;
   TempCellArray{1,15}=ZscoreDFFCornerZoneSigMean;
   TempCellArray{1,16}=DFFCornerZoneSigCounts;
   TempCellArray{1,17}=CornerZoneSignalRate;
   TempCellArray{1,18}=DFFCornerZoneStdev;
   TempCellArray{1,19}=CornerZoneCounts;
  
   TempCellArray{1,20}=DFFNonZoneMean;
   TempCellArray{1,21}=DFFNonZoneSigMean;
   TempCellArray{1,22}=ZscoreDFFNonZoneSigMean;
   TempCellArray{1,23}=DFFNonZoneSigCounts;
   TempCellArray{1,24}=NonZoneSignalRate;
   TempCellArray{1,25}=DFFNonZoneStdev;
   TempCellArray{1,26}=NonZoneCounts;
   
 
   TempCellArray{1,27}=ZonePeakHeight;
   TempCellArray{1,28}=ZscoreZonePeakHeight;
   TempCellArray{1,29}=RateZonePeakCounts;
   TempCellArray{1,30}=ZonePeakCounts;
  
   
   TempCellArray{1,31}=CornerPeakHeight;
   TempCellArray{1,32}=ZscoreCornerZonePeakHeight;
   TempCellArray{1,33}=RateCornerPeakCounts;
   TempCellArray{1,34}=CornerPeakCounts;
   
   TempCellArray{1,35}=NonZonePeakHeight;
   TempCellArray{1,36}=ZscoreNonZonePeakHeight;
   TempCellArray{1,37}=RateNonZonePeakCounts;
   TempCellArray{1,38}=NonZonePeakCounts;
   TempCellArray{1,39}=reply1;

   

   
    
%    TempCellLabel={'basename1' 'basename' 'NSegment' 'datastart' 'dataend' 'DFFZoneMean' 'DFFZoneStdev' 'DFFZoneCounts' 'DFFCornerZoneMean' 'DFFCornerZoneStdev' 'DFFCornerZoneCounts' 'DFFNonZoneMean' 'DFFNonZoneStdev' 'DFFNonZoneCounts' 'ZonePeakCounts' 'NormalizedZonePeakCounts' 'ZonePeakHeight' 'CornerZonePeakCounts' 'NormalizedCornerZonePeakCounts' 'CornerZonePeakHeight' 'NonZonePeakCounts' 'NormalizedNonZonePeakCounts' 'NonZonePeakHeight'};
   TempCellLabel={"basename1" "basename" "NSegment" "datastart" "dataend" "DFFZoneMean" "DFFZoneSigMean" "ZscoreDFFZoneSigMean" "DFFZoneSigCounts" "ZoneSignalRate" "DFFZoneStdev" "ZoneCounts" "DFFCornerZoneMean" "DFFCornerZoneSigMean" "ZscoreDFFCornerZoneSigMean" "DFFCornerZoneSigCounts"  "CornerZoneSignalRate" "DFFCornerZoneStdev" "CornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneSigMean" "ZscoreDFFNonZoneSigMean" "DFFNonZoneSigCounts" "NonZoneSignalRate" "DFFNonZoneStdev" "NonZoneCounts" "ZonePeakHeight" "ZscoreZonePeakHeight" "RateZonePeakCounts" "ZonePeakCounts" "CornerPeakHeight" "ZscoreCornerZonePeakHeight" "RateCornerPeakCounts" "CornerPeakCounts" "NonZonePeakHeight" "ZscoreNonZonePeakHeight" "RateNonZonePeakCounts" "NonZonePeakCounts" "TypeOfAnalysisControl" };
  
   TempCellLabel1=["basename1" "basename" "NSegment" "datastart" "dataend" "DFFZoneMean" "DFFZoneSigMean" "ZscoreDFFZoneSigMean" "DFFZoneSigCounts" "ZoneSignalRate" "DFFZoneStdev" "ZoneCounts" "DFFCornerZoneMean" "DFFCornerZoneSigMean" "ZscoreDFFCornerZoneSigMean" "DFFCornerZoneSigCounts"  "CornerZoneSignalRate" "DFFCornerZoneStdev" "CornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneSigMean" "ZscoreDFFNonZoneSigMean" "DFFNonZoneSigCounts" "NonZoneSignalRate" "DFFNonZoneStdev" "NonZoneCounts" "ZonePeakHeight" "ZscoreZonePeakHeight" "RateZonePeakCounts" "ZonePeakCounts" "CornerPeakHeight" "ZscoreCornerZonePeakHeight" "RateCornerPeakCounts" "CornerPeakCounts" "NonZonePeakHeight" "ZscoreNonZonePeakHeight" "RateNonZonePeakCounts" "NonZonePeakCounts" "TypeOfAnalysisControl"];

%    datafilename=['FiberDataExport4' '.csv'];
%     if ~isempty(datafilename);
%         xlswrite(datafilename,TempCellLabel,'Sheet1');
%     else
%     end
% 
%     [success,message]=xlsappend(datafilename,TempCellArray);
   TableTempCellLabel=cell2table(TempCellLabel, "VariableNames", TempCellLabel1);
   datafilename=['125_xs-SLR_TEST_summary_all_080525' '.csv'];
    if isfile(datafilename);
      
    else
    writetable(TableTempCellLabel,datafilename);
    end
    
     CellTable=table2cell(readtable(datafilename));
     CellTableNew=[CellTable;TempCellArray];
     TableNew=cell2table(CellTableNew, "VariableNames", TempCellLabel1);
     writetable(TableNew,datafilename);

   


    

     
%%  
%reply = input('Do you want PSTH? y/n [y]:','s');
   
   %if isempty(reply);
     %reply = 'n';
   %end
   
% Creating a PSTH with Df/f and zscore ( 1s in non zone and 2s in zone 3/14/25 PKN % Updated 3/17/25 
% Assume data is stored in an array called data
% Columns are as follows:
% Column 1 - Frame Number
% Column 2 - df/f
% Column 3 - Z-score
% Column 4 - Zone indicator (0 or 1)
% Column 5 - (Unused, optional)
% Column 6 - Zone indicator (0 or 1)

Frames=zeros(lengthFiberData, 4);
Frames(:,1)= combinedFibVidData(:,2);
Frames(:,2)= combinedFibVidData(:,6);
Frames(:,3)= combinedFibVidData(:,7);
Frames(:,4:6)= combinedFibVidData(:,8:10);


%%INTERACTION ZONE 
numFrames = size(Frames, 1);  % Get the total number of frames
FilteredFrames = {};        % Initialize an empty cell array to store valid event sets
i = 1;  % Start from frame 1

while i <= numFrames
    % Find the next occurrence of '1' in column 4 after the current frame
    if Frames(i, 4) == 1
        % Check if there are at least 10 preceding frames where column 6 is 1
        if i > 10
            prevFrames = Frames(i-10:i-1, 6);  % Extract the 10 preceding frames from column 6
            
            % Check if all 10 preceding frames have 1 in column 6
            if all(prevFrames == 1)
                % Initialize event frames with the 10 preceding frames
                eventFrames = Frames(i-10:i, :);

                % Now add the succeeding frames where column 4 is 1 (up to 20 frames)
                j = i;
                count = 0;
                while j <= numFrames && Frames(j, 4) == 1 && count < 20
                    eventFrames = [eventFrames; Frames(j, :)];  % Add the frame to the event
                    j = j + 1;
                    count = count + 1;
                end

                % Ensure the event has exactly 30 frames (10 preceding + 20 succeeding)
                % If there are more than 20 succeeding frames, we limit the event to the first 20 frames
                if size(eventFrames, 1) < 30
                    % Add the event frames as a separate event set to FilteredFrames
                    FilteredFrames{end+1} = eventFrames;
                else
                    % If eventFrames is larger than 30 frames (i.e., more than 20 succeeding frames)
                    % we just take the first 30 frames
                    FilteredFrames{end+1} = eventFrames(1:30, :);
                end

                % Skip to the next occurrence of 1 (next event) beyond the current event's end
                i = j + 1;
            else
                % If the preceding frames condition is not met, skip this event
                i = i + 1;
            end
        else
            % Skip if not enough preceding frames (less than 10)
            i = i + 1;
        end
    else
        % Move to the next frame if no '1' is found in column 4
        i = i + 1;
    end
end

% Display the total number of valid event sets
disp(['Total number of valid event sets: ', num2str(length(FilteredFrames))]);

% Display the filtered frames for the first event (if any)
if ~isempty(FilteredFrames)
    disp('Filtered event set (first one):');
    disp(FilteredFrames{1});
end




% Initialize an array to hold only valid events with at least 30 frames
validFilteredFrames = {};

% Loop through each event in the FilteredFrames cell array
for idx = 1:length(FilteredFrames)
    eventData = FilteredFrames{idx};
    
    % If the event has fewer than 30 frames, skip this event
    if size(eventData, 1) < 30
        continue;  % Skip this event
    end
    
    % Store valid event
    validFilteredFrames{end+1} = eventData;
end

% Update FilteredFrames to only contain valid events
FilteredFrames = validFilteredFrames;

% Create a figure for plotting with 2 subplots
%H5=figure('Name',basename2);
%figure(H5);

% Initialize matrices to store Z-scores and df/F of valid events
allZData = [];
allDF_FData = [];
maxFrames = 30; % Fixed to 30 since only valid events are included

% -------------------- First Plot: Z-score --------------------
%subplot(2,2,1);  % Create the first subplot (Z-score)
%title('Z-score for Interaction Zone');
%xlabel('Time (s)');
%ylabel('Z-score');
%hold on;

% Loop through each event in the valid FilteredFrames cell array
for idx = 1:length(FilteredFrames)
    eventData = FilteredFrames{idx};
    
    % Extract Z-score (column 3)
    zScore = eventData(:, 3);
    
    % Convert frame indices to time
    timeAxis = (1:30) - 10;  % Shift frames so that Frame 10 is at time 0
    timeAxis = timeAxis / 10;  % Convert to seconds
    
    % Plot individual Z-score events (gray)
    %plot(timeAxis, zScore, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1); 

    % Store event Z-scores for averaging
    allZData = [allZData; zScore'];
end

% Compute the mean Z-score
meanZScore = mean(allZData, 1);

% Plot the averaged Z-score as a bold black line
%plot(timeAxis, meanZScore, 'k', 'LineWidth', 2);
%yLimits = ylim; % Get current y-axis limits
    %plot([0 0], [yLimits(1), yLimits(2)], 'r', 'LineWidth', 1);
    
    % Re-fetch limits in case they changed
    %yLimits = ylim; 
    
   %Adding text
    %text(0, yLimits(2) * 0.98, 'Zone entry', 'Color', 'r', 'FontSize', 6, ...
        %'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        %'VerticalAlignment', 'bottom');
    %writematrix(meanZScore, '5611(0G)_TEST_meanZscore_Interaction Zone_032625.txt')

%hold off;

% -------------------- Second Plot: df/F --------------------
%subplot(2,2,2);  % Create the second subplot (df/F)
%title('df/F for Interaction Zone');
%xlabel('Time (s)');
%ylabel('df/F');
%hold on;

% Loop through each event again for df/F
for idx = 1:length(FilteredFrames)
    eventData = FilteredFrames{idx};
    
    % Extract df/F (column 2)
    df_F = eventData(:, 2);
    
    % Plot individual df/F events (light blue)
    plot(timeAxis, df_F, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1); 

    % Store event df/F values for averaging
    allDF_FData = [allDF_FData; df_F'];
end

% Compute the mean df/F
meanDF_F = mean(allDF_FData, 1);

% Plot the averaged df/F as a bold black line
%plot(timeAxis, meanDF_F, 'b', 'LineWidth', 2); 
%yLimits = ylim; % Get current y-axis limits
    %plot([0 0], [yLimits(1), yLimits(2)], 'r', 'LineWidth', 1);
    
    % Re-fetch limits in case they changed
    %yLimits = ylim; 
    
    % Add text annotation slightly above the top of the plot
    %text(0, yLimits(2) * 0.98, 'Zone entry', 'Color', 'r', 'FontSize', 6, ...
        %'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        %'VerticalAlignment', 'bottom');
 %writematrix(meanDF_F, '5611(0G)_TEST_meanDF_F_Interaction Zone_032625.txt')

%%CORNER ZONE 

numFrames = size(Frames, 1);  % Get the total number of frames
FilteredFramesCZ = {};        % Initialize an empty cell array to store valid event sets
i = 1;  % Start from frame 1
while i <= numFrames
    % Find the next occurrence of '1' in column 5 after the current frame
    if Frames(i, 5) == 1
        % Check if there are at least 10 preceding frames where column 6 is 1
        if i > 10
            prevFramesCZ = Frames(i-10:i-1, 6);  % Extract the 10 preceding frames from column 6
            
            % Check if all 10 preceding frames have 1 in column 6
            if all(prevFramesCZ == 1)
                % Initialize event frames with the 10 preceding frames
                eventFramesCZ = Frames(i-10:i, :);

                % Now add the succeeding frames where column 5 is 1 (up to 10 frames)
                j = i;
                count = 0;
                while j <= numFrames && Frames(j, 5) == 1 && count < 10
                    eventFramesCZ = [eventFramesCZ; Frames(j, :)];  % Add the frame to the event
                    j = j + 1;
                    count = count + 1;
                end

                % Ensure the event has exactly 20 frames (10 preceding + 10 succeeding)
                % If there are more than 20 succeeding frames, we limit the event to the first 20 frames
                if size(eventFramesCZ, 1) < 20
                    % Add the event frames as a separate event set to FilteredFrames
                    FilteredFramesCZ{end+1} = eventFramesCZ;
                else
                    % If eventFrames is larger than 20 frames (i.e., more than 10 succeeding frames)
                    % we just take the first 20 frames
                    FilteredFramesCZ{end+1} = eventFramesCZ(1:20, :);
                end

                % Skip to the next occurrence of 1 (next event) beyond the current event's end
                i = j + 1;
            else
                % If the preceding frames condition is not met, skip this event
                i = i + 1;
            end
        else
            % Skip if not enough preceding frames (less than 10)
            i = i + 1;
        end
    else
        % Move to the next frame if no '1' is found in column 4
        i = i + 1;
    end
end

% Display the total number of valid event sets
disp(['Total number of valid event sets: ', num2str(length(FilteredFramesCZ))]);

% Display the filtered frames for the first event (if any)
if ~isempty(FilteredFramesCZ)
    disp('Filtered event set (first one):');
    disp(FilteredFramesCZ{1});
end


% Initialize an array to hold only valid events with at least 30 frames
validFilteredFramesCZ = {};

% Loop through each event in the FilteredFrames cell array
for idx = 1:length(FilteredFramesCZ)
    eventDataCZ = FilteredFramesCZ{idx};
    
    % If the event has fewer than 20 frames, skip this event
    if size(eventDataCZ, 1) < 20
        continue;  % Skip this event
    end
    
    % Store valid event
    validFilteredFramesCZ{end+1} = eventDataCZ;
end

% Update FilteredFrames to only contain valid events
FilteredFramesCZ = validFilteredFramesCZ;

% Initialize matrices to store Z-scores and df/F of valid events
allZDataCZ = [];
allDF_FDataCZ = [];
maxFrames = 20; % Fixed to 20 since only valid events are included

% -------------------- First Plot: Z-score --------------------
%subplot(2,2,3);  % Create the first subplot (Z-score)
%title('Z-score for CORNERZone');
%xlabel('Time (s)');
%ylabel('Z-score');
%hold on;

% Loop through each event in the valid FilteredFrames cell array
for idx = 1:length(FilteredFramesCZ)
    eventDataCZ = FilteredFramesCZ{idx};
    
    % Extract Z-score (column 3)
    zScoreCZ = eventDataCZ(:, 3);
    
    % Convert frame indices to time
    timeAxis = (1:20) - 10;  % Shift frames so that Frame 10 is at time 0
    timeAxis = timeAxis / 10;  % Convert to seconds
    
    % Plot individual Z-score events (gray)
    plot(timeAxis, zScoreCZ, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1); 

    % Store event Z-scores for averaging
    allZDataCZ = [allZDataCZ; zScoreCZ'];
end
if isempty(allZDataCZ);
    warning('No valid events found. Skipping mean calculation and plotting');
else

    
    %Compute the mean Z-score
meanZScoreCZ = mean(allZDataCZ, 1);

% Plot the averaged Z-score as a bold black line
%plot(timeAxis, meanZScoreCZ, 'k', 'LineWidth', 2); 
%yLimits = ylim; % Get current y-axis limits
    %plot([0 0], [yLimits(1), yLimits(2)], 'r', 'LineWidth', 1);
    
    % Re-fetch limits in case they changed
    %yLimits = ylim; 
    
    % Add text annotation slightly above the top of the plot
    %text(0, yLimits(2) * 0.98, 'Zone entry', 'Color', 'r', 'FontSize', 6, ...
        %'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        %'VerticalAlignment', 'bottom');

%hold off;

% -------------------- Second Plot: df/F --------------------
%subplot(2,2,4);  % Create the second subplot (df/F)
%title('df/F for CORNERZone');
%xlabel('Time (s)');
%ylabel('df/F');
%hold on;

% Loop through each event again for df/F
for idx = 1:length(FilteredFramesCZ)
    eventDataCZ = FilteredFramesCZ{idx};
    
    % Extract df/F (column 2)
    df_FCZ = eventDataCZ(:, 2); 
    
    % Plot individual df/F events (light blue)
    %plot(timeAxis, df_FCZ, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1); 

    % Store event df/F values for averaging
    %allDF_FDataCZ = [allDF_FDataCZ; df_F'];
end

% Compute the mean df/F
meanDF_FCZ = mean(allDF_FDataCZ, 1);

% Plot the averaged df/F as a bold black line
%plot(timeAxis, meanDF_FCZ, 'k', 'LineWidth', 2); 

%yLimits = ylim; % Get current y-axis limits
   %plot([0 0], [yLimits(1), yLimits(2)], 'r', 'LineWidth', 1);
    
    % Re-fetch limits in case they changedSigtable
    %yLimits = ylim; 
    
    % Add text annotation slightly above the top of the plot
    %text(0, yLimits(2) * 0.98, 'Zone entry', 'Color', 'r', 'FontSize', 6, ...
        %'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        %'VerticalAlignment', 'bottom');
%hold off;
%figurefilename6=[basename1,'PSTH.png'];
   %saveas(H5,figurefilename6);

writematrix(allDF_FData, '125_xs-SLR_TEST_Novel_DFF_PSTH.csv'); 
writematrix(allDF_FDataCZ, '125_xs-SLR_TEST_Familiar_DFF_PSTH.csv'); 
writematrix(allZData, '125_xs-SLR_TEST_Novel_Zscore_PSTH.csv'); 
writematrix(allZDataCZ, '125_xs-SLR_TEST_Familiar_Zscore_PSTH.csv'); 
   
   

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
   partialAnalysisData(:,7:9)=finalCombinedData(startpoint:endpoint,5:7);
   
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
    subplot(5,2,1)
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
    
    subplot(6,2,[3,4])
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
    subplot(5,2,2)
    plot(sigmaplotF);
    
    
    savefilename0=[basename1 '_' num2str(startpoint) '_' num2str(endpoint) '-ALL.txt'];
    dlmwrite(savefilename0,partialAnalysisData);
 
    MAD=median(abs(DeltaFoverFF-meanFF));
    [pksPA,locsPA]=findpeaks(DeltaFoverFF-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);
    [pks02,locs02]=findpeaks(DeltaFoverFF-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0);%defining "no thresold" % 08/04/25 PKN
    
    %[pks2,locs2]=findpeaks(sigmaplotF,'MinPeakHeight',2,'MinPeakDistance',3);
    numpks2=length(locsPA);
    numpks02=length(locs02);% 08/04/25 PKN
    
       
    peaktable2=zeros(numpks2,3);%initialize peaktable 
    peaktable2(:,1)=locsPA;
    peaktable2(:,2)=pksPA;
    
    dff0pks2=DeltaFoverFF(locsPA);
    peaktable2(:,3)=dff0pks2;
    peaktable2(:,4) = fiberZScores(locsPA,2);

    %initialize a table containing all signals before threshold % 08/04/25 PKN
    PASigtable=zeros(numpks02,3);
    PASigtable(:,1)=locs02;
    PASigtable(:,2)=pks02;
    dff0pks02=DeltaFoverF(locs02);
    PASigtable(:,3)=dff0pks02;
    PASigtable(:,4) = fiberZScores(locs02,2);
    PASigtable(:,5) = combinedFibVidData(locs02,3);
    locs= combinedFibVidData(locs0,3); % extracts corresponding time for an occurring signal
    
    
    figure(H1);
    subplot(5,2,[5,6])
    plot(partialAnalysisData(:,6),'Color','black');
    hold on 
    for ii=1:numpks2;
        scatter(locsPA(ii),dff0pks2(ii),5,'blue');
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
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.035*CCCC(:,2),'.','Color',[0.145, 0.537, 0.741]);
    ylim([-0.050, max(ylim)+0.005]);
    
    %{
%second graph with time in s as x-axis (time from FP) %08/04/2025 PKN
    figure(H1);
    subplot(5,2,[7,8])
    plot(partialAnalysisData(:,3), partialAnalysisData(:,6),'Color','black');
    hold on 
    for ii=1:numpks2;
        scatter(locsPA(ii),dff0pks2(ii),5,'blue');
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
        counter=counter+1;
        else
        end
    end
    plot(CCCC(:,1),-0.035*CCCC(:,2),'.','Color',[0.145, 0.537, 0.741]);
    ylim([-0.050, max(ylim)+0.005]);
    %}
    
    
    
    figurefilename2=[basename1 '_' num2str(startpoint) '_' num2str(endpoint) '-.png'];
    saveas(H1,figurefilename2);

    
    locs3=[];
    %creating separate events for averaging (CHANGE ALL OF THE VARIABLES)

   finalCombinedDataTempPA(:,1) = bwlabel(logical(partialAnalysisData(:,7)));
   maxNumEvents1 = max(finalCombinedDataTempPA(:,1));
   

    if maxNumEvents1 > 0
    TempCellArrayEventPA1=[];    
    TempCellArrayEventPA1{maxNumEvents1,10}=[];     
        for i = 1: maxNumEvents1
        locs3 = find(finalCombinedDataTempPA(:,1) == i);
        if numel(locs3) > 5;
            FiberDataArrayPA.Active{i} = partialAnalysisData(locs3, 4);
            finalCombinedDataTempPA1 = partialAnalysisData(locs3, 4);
            [pksPAEvent1,locsPAEvent1]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);%May23,2023
            TempCellArrayEventPA1{i,1}=basename1;
            TempCellArrayEventPA1{i,2}=basename2;
            TempCellArrayEventPA1{i,3}="Active";
            TempCellArrayEventPA1{i,4}=mean(partialAnalysisData(locs3, 6));%May23,2023
            TempCellArrayEventPA1{i,5}=mean(sigmaplotF(locs3));% from line 802  % mean Z-score 
            TempCellArrayEventPA1{i,6}=numel(locsPAEvent1);
            TempCellArrayEventPA1{i,7}=mean(finalCombinedDataTempPA1(locsPAEvent1));
            TempCellArrayEventPA1{i,8} = str2num(reply1);
            TempCellArrayEventPA1{i,9}=numel(finalCombinedDataTempPA1);
            sigmaplotFtemp1=sigmaplotF(locs3); %HT 02/02/2024
            TempCellArrayEventPA1{i,10}=mean(sigmaplotFtemp1(locsPAEvent1)); % Mean Peak Height Z-Score
            TempCellArrayEventPA1{i,11}= numel(sigmaplotFtemp1(locsPAEvent1));

            [pksPAEvent01,locsPAEvent01]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %02/27/2025 PKN
            TempCellArrayEventPA1{i,12}= numel(locsPAEvent01);%02/27/2025 PKN
            TempCellArrayEventPA1{i,13}= mean(finalCombinedDataTempPA1(locsPAEvent01));%03/03/2025 PKN
            TempCellArrayEventPA1{i,14}=mean(sigmaplotFtemp1(locsPAEvent01));%03/03/2025 PKN
    
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
            TempCellArrayEventPA2{i,3}="Corner";
            TempCellArrayEventPA2{i,4}=mean(partialAnalysisData(locs3, 6)); %May 2023
            TempCellArrayEventPA2{i,5}=mean(sigmaplotF(locs3));% from line 802  % mean Z-score 
            TempCellArrayEventPA2{i,6}=numel(locsPAEvent2);
            TempCellArrayEventPA2{i,7}=mean(finalCombinedDataTempPA2(locsPAEvent2));
            TempCellArrayEventPA2{i,8} = str2num(reply1);
            TempCellArrayEventPA2{i,9}=numel(finalCombinedDataTempPA2);
            sigmaplotFtemp2=sigmaplotF(locs3); %HT 02/02/2024
            TempCellArrayEventPA2{i,10}=mean(sigmaplotFtemp2(locsPAEvent2)); % Mean Peak Height Z-Score
            TempCellArrayEventPA2{i,11}= numel(sigmaplotFtemp2(locsPAEvent2));
    
            [pksPAEvent02,locsPAEvent02]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %02/27/2025 PKN
            TempCellArrayEventPA2{i,12}= numel(locsPAEvent02);%02/27/2025 PKN
            TempCellArrayEventPA2{i,13}= mean(finalCombinedDataTempPA2(locsPAEvent02));%03/03/2025 PKN
            TempCellArrayEventPA2{i,14}=mean(sigmaplotFtemp2(locsPAEvent02));%03/03/2025 PKN

        end
        end
    else 
        
    TempCellArrayEventPA2=[];     
    end

    finalCombinedDataTempPA(:,3) = bwlabel(logical(partialAnalysisData(:,9)));
    maxNumEvents3 = max(finalCombinedDataTempPA(:,3));
    

    if maxNumEvents3 > 0
    TempCellArrayEventPA3=[];      
    TempCellArrayEventPA3{maxNumEvents3,10}=[];     
        for i = 1: maxNumEvents3
        locs3 = find(finalCombinedDataTempPA(:,3) == i);
         if numel(locs3) > 5
            FiberDataArrayPA.NonActive{i} = partialAnalysisData(locs3, 4);
            finalCombinedDataTempPA3 = partialAnalysisData(locs3, 4);
            [pksPAEvent3,locsPAEvent3]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',2*MAD,'MinPeakDistance',3);%May23,2023
            TempCellArrayEventPA3{i,1}=basename1;
            TempCellArrayEventPA3{i,2}=basename2;
            TempCellArrayEventPA3{i,3}="NonActive";
            TempCellArrayEventPA3{i,4}=mean(partialAnalysisData(locs3, 6)); %May 2023
            TempCellArrayEventPA3{i,5}=mean(sigmaplotF(locs3));% from line 802  % mean Z-score 
            TempCellArrayEventPA3{i,6}=numel(locsPAEvent3);
            TempCellArrayEventPA3{i,7}=mean(finalCombinedDataTempPA3(locsPAEvent3));
            TempCellArrayEventPA3{i,8} = str2num(reply1);
            TempCellArrayEventPA3{i,9}=numel(finalCombinedDataTempPA3);
            sigmaplotFtemp3=sigmaplotF(locs3); %HT 02/02/2024
            TempCellArrayEventPA3{i,10}=mean(sigmaplotFtemp3(locsPAEvent3)); % Mean Peak Height Z-Score
            TempCellArrayEventPA3{i,11}= numel(sigmaplotFtemp3(locsPAEvent3));

            [pksPAEvent03,locsPAEvent03]=findpeaks(partialAnalysisData(locs3, 6)-meanFF,'MinPeakHeight',0.1*MAD,'MinPeakDistance',0); %02/27/2025 PKN
            TempCellArrayEventPA3{i,12}= numel(locsPAEvent03);%02/27/2025 PKN
            TempCellArrayEventPA3{i,13}= mean(finalCombinedDataTempPA3(locsPAEvent03));%03/03/2025 PKN
            TempCellArrayEventPA3{i,14}=mean(sigmaplotFtemp3(locsPAEvent03));%03/03/2025 PKN

         end
        end
    else 
    TempCellArrayEventPA3=[];     
    end

   TempCellArrayAllEventsPA=[];
   TempCellArrayAllEventsPA = [TempCellArrayEventPA1 ; TempCellArrayEventPA2; TempCellArrayEventPA3];

   TempCellLabelEventsPA={"basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z", "numofSig(zs)" "numSig" "meanSigDFF" "MeanSigDFF in Z"};
   TempCellLabelEventsPA1=["basename1" "basename" "EventType" "MeanDFF" "MeanZScore" "numPeaks" "meanPeakHeight" "TypeOfAnalysisControl" "duration" "meanPeakHeight in Z" "numofSig(zs)" "numSig" "meanSigDFF" "MeanSigDFF in Z"];


   TableTempCellLabelEventsPA=cell2table(TempCellLabelEventsPA, "VariableNames", TempCellLabelEventsPA1);
   datafilenameEventsPA=['125_xs-SLR_TEST_events_partial_080425' '.csv'];
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
    
    TotalActiveCountF=sum(partialAnalysisData(:,7));
    TempArrayF=zeros(TotalActiveCountF,1);
    
    TotalCornerCountF=sum(partialAnalysisData(:,8));
    TempArray2F=zeros(TotalCornerCountF,1);
    
    TotalNonActiveCountF=sum(partialAnalysisData(:,9));
    TempArray3F=zeros(TotalNonActiveCountF,1);
    
    nn=1; %counter 
    nnn=1; %counter 
    nnnn=1; %counter 
    for iii=1:lengthF 
     if partialAnalysisData(iii,7)==1;
         TempArrayF(nn,1)=partialAnalysisData(iii,6);
         nn=nn+1;
     elseif partialAnalysisData(iii,8)==1
         TempArray2F(nnn,1)=partialAnalysisData(iii,6);
         nnn=nnn+1;
     else
         TempArray3F(nnnn,1)=partialAnalysisData(iii,6);
         nnnn=nnnn+1;
     end
    end
    
    %Finding DF/F signal amplitude (before thresold) %08/04/25 PKN 
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks02;
        NNN=PASigtable(iii,1);
        if finalCombinedData(NNN,5)==1;
            PAActiveZoneSig(nn,1)=PASigtable(iii,3);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            PACornerZoneSig(nnn,1)=PASigtable(iii,3);
            nnn=nnn+1;
        else
            PANonActiveZoneSig(nnnn,1)=PASigtable(iii,3);
            nnnn=nnnn+1;
        end 
    end

%Finding Z-score of DF/F signal amplitude (before thresold) %08/04/25 PKN 
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks02;
        NNN=Sigtable(iii,1);
        if finalCombinedData(NNN,5)==1;
            PAZscoreActiveZoneSig(nn,1)=PASigtable(iii,4);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            PAZscoreCornerZoneSig(nnn,1)=PASigtable(iii,4);
            nnn=nnn+1;
        else
            PAZscoreNonActiveZoneSig(nnnn,1)=PASigtable(iii,4);
            nnnn=nnnn+1;
        end 
    end

    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;
    for iii=1:numpks2;
        NNN=peaktable2(iii,1);
        if partialAnalysisData(NNN,7)==1;
            PAActiveZonePeakF(nn,1)=peaktable2(iii,3);
            nn=nn+1;
        elseif partialAnalysisData(NNN,8)==1;
            PACornerZonePeakF(nnn,1)=peaktable2(iii,3);
            nnn=nnn+1;
        else
            PANonActiveZonePeakF(nnnn,1)=peaktable2(iii,3);
            nnnn=nnnn+1;
        end 
    end
    % Finding Z-score of peak amplitude  %08/04/25 PKN     
    nn=1;% counter 
    nnn=1; % counter 
    nnnn=1;

    for iii=1:numpks2;
        NNN=peaktable(iii,1);
        if finalCombinedData(NNN,5)==1;
            PAZscoreActiveZonePeak(nn,1)=peaktable2(iii,4);
            nn=nn+1;
        elseif finalCombinedData(NNN,6)==1;
            PAZscoreCornerZonePeak(nnn,1)=peaktable2(iii,4);
            nnn=nnn+1;
        else
            PAZscoreNonActiveZonePeak(nnnn,1)=peaktable2(iii,4);
            nnnn=nnnn+1;
        end 
    end  
   
   if nn==1;
       PAActiveZonePeakF=[];
   end
   
   if nnn==1;
       PACornerZonePeakF=[];
   else
   end
   
   if nnnn==1;
       PANonActiveZonePeakF=[];
   else
   end

   PACornerZoneSig = []; 

   if ~isempty(PACornerZoneSig)
    % The array has data, so calculate the mean.
    DFFCornerZoneSigMean = mean(PACornerZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
    DFFCornerZoneSigMean = NaN; 
   end

   PAZscoreActiveZoneSig = []; 

   if ~isempty(PAActiveZoneSig)
    % The array has data, so calculate the mean.
    ZscoreDFFZoneSigMean = mean(PAZscoreActiveZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
    ZscoreDFFZoneSigMean = NaN; 
   end

   PAZscoreCornerZoneSig = []; 

   if ~isempty(PACornerZoneSig)
    % The array has data, so calculate the mean.
    ZscoreDFFCornerZoneSigMean = mean(PAZscoreCornerZoneSig);
    else
    % The array is empty. Set the mean to a default value.
    % Using NaN (Not a Number) is often a good practice,
    % but you could also use 0 if that makes sense for your data.
    ZscoreDFFCornerZoneSigMean = NaN; 
   end
    
   DFFZoneMean=mean(TempArrayF);
   DFFZoneSigMean=mean(PAActiveZoneSig);
   ZscoreDFFZoneSigMean=mean(PAZscoreActiveZoneSig);
   DFFZoneStdev=std(TempArrayF);
   ZoneCounts=size(TempArrayF,1);
   DFFZoneSigCounts=size(PAActiveZoneSig,1);
   ZoneSignalRate=DFFZoneSigCounts/ZoneCounts*20;% peak#/sec

   DFFCornerZoneMean=mean(TempArray2F);
   DFFCornerZoneSigMean=mean(PACornerZoneSig);
   ZscoreDFFCornerZoneSigMean= mean(PAZscoreCornerZoneSig);
   DFFCornerZoneStdev=std(TempArray2F);
   CornerZoneCounts=size(TempArray2F,1); 
   DFFCornerZoneSigCounts=size(PACornerZoneSig,1);
   CornerZoneSignalRate=DFFCornerZoneSigCounts/CornerZoneCounts*20;% peak#/sec 

   DFFNonZoneMean=mean(TempArray3F);
   DFFNonZoneSigMean=mean(PANonActiveZoneSig);
   ZscoreDFFNonZoneSigMean=mean(ZscoreNonActiveZoneSig);
   DFFNonZoneStdev=std(TempArray3F);
   NonZoneCounts=size(TempArray3F,1); 
   DFFNonZoneSigCounts=size(PANonActiveZoneSig,1);
   NonZoneSignalRate=DFFNonZoneSigCounts/NonZoneCounts*20;% peak#/sec

   ZonePeakCounts=size(PAActiveZonePeakF,1);
   RateZonePeakCounts=ZonePeakCounts/ZoneCounts*20;% peak#/sec
   ZonePeakHeight=mean(PAActiveZonePeakF);
   ZscoreZonePeakHeight=mean(PAZscoreActiveZonePeak);

   CornerPeakCounts=size(PACornerZonePeakF,1);
   RateCornerPeakCounts=CornerPeakCounts/CornerZoneCounts*20;
   CornerPeakHeight=mean(PACornerZonePeakF);
   ZscoreCornerZonePeakHeight=mean(PAZscoreCornerZonePeak); 



   NonZonePeakCounts=size(PANonActiveZonePeakF,1);
   RateNonZonePeakCounts=NonZonePeakCounts/NonZoneCounts*20;
   NonZonePeakHeight=mean(PANonActiveZonePeakF);
   ZscoreNonZonePeakHeight=mean(PAZscoreNonActiveZonePeak); 

    

    
   TempCellArray{1,39}=[]; 
   TempCellArray{1,1}=basename1;
   TempCellArray{1,2}=basename2;
   
   TempCellArray{1,3}=NSegment;
   
   TempCellArray{1,4}=startpoint;
   TempCellArray{1,5}=endpoint;
  
   TempCellArray{1,6}=DFFZoneMean;
   TempCellArray{1,7}=DFFZoneSigMean;
   TempCellArray{1,8}=ZscoreDFFZoneSigMean;
   TempCellArray{1,9}=DFFZoneSigCounts;
   TempCellArray{1,10}=ZoneSignalRate;
   TempCellArray{1,11}=DFFZoneStdev;
   TempCellArray{1,12}=ZoneCounts;
   
   TempCellArray{1,13}=DFFCornerZoneMean;
   TempCellArray{1,14}=DFFCornerZoneSigMean;
   TempCellArray{1,15}=ZscoreDFFCornerZoneSigMean;
   TempCellArray{1,16}=DFFCornerZoneSigCounts;
   TempCellArray{1,17}=CornerZoneSignalRate;
   TempCellArray{1,18}=DFFCornerZoneStdev;
   TempCellArray{1,19}=CornerZoneCounts;
  
   TempCellArray{1,20}=DFFNonZoneMean;
   TempCellArray{1,21}=DFFNonZoneSigMean;
   TempCellArray{1,22}=ZscoreDFFNonZoneSigMean;
   TempCellArray{1,23}=DFFNonZoneSigCounts;
   TempCellArray{1,24}=NonZoneSignalRate;
   TempCellArray{1,25}=DFFNonZoneStdev;
   TempCellArray{1,26}=NonZoneCounts;
   
 
   TempCellArray{1,27}=ZonePeakHeight;
   TempCellArray{1,28}=ZscoreZonePeakHeight;
   TempCellArray{1,29}=RateZonePeakCounts;
   TempCellArray{1,30}=ZonePeakCounts;
  
   
   TempCellArray{1,31}=CornerPeakHeight;
   TempCellArray{1,32}=ZscoreCornerZonePeakHeight;
   TempCellArray{1,33}=RateCornerPeakCounts;
   TempCellArray{1,34}=CornerPeakCounts;
   
   TempCellArray{1,35}=NonZonePeakHeight;
   TempCellArray{1,36}=ZscoreNonZonePeakHeight;
   TempCellArray{1,37}=RateNonZonePeakCounts;
   TempCellArray{1,38}=NonZonePeakCounts;
   TempCellArray{1,39}=reply1;

    
%    TempCellLabel={'basename1' 'basename' 'NSegment' 'datastart' 'dataend' 'DFFZoneMean' 'DFFZoneStdev' 'DFFZoneCounts' 'DFFCornerZoneMean' 'DFFCornerZoneStdev' 'DFFCornerZoneCounts' 'DFFNonZoneMean' 'DFFNonZoneStdev' 'DFFNonZoneCounts' 'ZonePeakCounts' 'NormalizedZonePeakCounts' 'ZonePeakHeight' 'CornerZonePeakCounts' 'NormalizedCornerZonePeakCounts' 'CornerZonePeakHeight' 'NonZonePeakCounts' 'NormalizedNonZonePeakCounts' 'NonZonePeakHeight'};
   TempCellLabel={"basename1" "basename" "NSegment" "datastart" "dataend" "DFFZoneMean" "DFFZoneSigMean" "ZscoreDFFZoneSigMean" "DFFZoneSigCounts" "ZoneSignalRate" "DFFZoneStdev" "ZoneCounts" "DFFCornerZoneMean" "DFFCornerZoneSigMean" "ZscoreDFFCornerZoneSigMean" "DFFCornerZoneSigCounts"  "CornerZoneSignalRate" "DFFCornerZoneStdev" "CornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneSigMean" "ZscoreDFFNonZoneSigMean" "DFFNonZoneSigCounts" "NonZoneSignalRate" "DFFNonZoneStdev" "NonZoneCounts" "ZonePeakHeight" "ZscoreZonePeakHeight" "RateZonePeakCounts" "ZonePeakCounts" "CornerPeakHeight" "ZscoreCornerZonePeakHeight" "RateCornerPeakCounts" "CornerPeakCounts" "NonZonePeakHeight" "ZscoreNonZonePeakHeight" "RateNonZonePeakCounts" "NonZonePeakCounts" "TypeOfAnalysisControl" };
  
   TempCellLabel1=["basename1" "basename" "NSegment" "datastart" "dataend" "DFFZoneMean" "DFFZoneSigMean" "ZscoreDFFZoneSigMean" "DFFZoneSigCounts" "ZoneSignalRate" "DFFZoneStdev" "ZoneCounts" "DFFCornerZoneMean" "DFFCornerZoneSigMean" "ZscoreDFFCornerZoneSigMean" "DFFCornerZoneSigCounts"  "CornerZoneSignalRate" "DFFCornerZoneStdev" "CornerZoneCounts" "DFFNonZoneMean" "DFFNonZoneSigMean" "ZscoreDFFNonZoneSigMean" "DFFNonZoneSigCounts" "NonZoneSignalRate" "DFFNonZoneStdev" "NonZoneCounts" "ZonePeakHeight" "ZscoreZonePeakHeight" "RateZonePeakCounts" "ZonePeakCounts" "CornerPeakHeight" "ZscoreCornerZonePeakHeight" "RateCornerPeakCounts" "CornerPeakCounts" "NonZonePeakHeight" "ZscoreNonZonePeakHeight" "RateNonZonePeakCounts" "NonZonePeakCounts" "TypeOfAnalysisControl"];

  %    datafilename=['FiberDataExport4' '.csv'];
%     if ~isempty(datafilename);
%         xlswrite(datafilename,TempCellLabel,'Sheet1');
%     else
%     end
% 
%     [success,message]=xlsappend(datafilename,TempCellArray);
   TableTempCellLabel=cell2table(TempCellLabel, "VariableNames", TempCellLabel1);
   datafilename=['125_xs-SLR_TEST_summary_PARTIAL_080525' '.csv'];
    if isfile(datafilename);
      
    else
    writetable(TableTempCellLabel,datafilename);
    end
    
     CellTable=table2cell(readtable(datafilename));
     CellTableNew=[CellTable;TempCellArray];
     TableNew=cell2table(CellTableNew, "VariableNames", TempCellLabel1);
     writetable(TableNew,datafilename);
 else
 end   
    
%%
   


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