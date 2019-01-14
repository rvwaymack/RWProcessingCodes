%% Burst calling based on slope threshold 

%Look at nc14 only
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);

for n=1:length(SpotDiff)
    %Limit to spots that exist for 3 or more frames
    if length(SpotDiff(n).Frame1) >=3 & length(SpotDiff(n).SpotOne)==length(ElapsedTime)
       %limit to 50minutes into nc14 (100 frames)
       if (length(ElapsedTime)-nc14)>=100
       SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:(nc14+100)),SpotDiff(n).SpotOne(nc14:(nc14+100)),0.1,'lowess');
       else
           ShortTrace='y';
           SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiff(n).SpotOne(nc14:end),0.1,'lowess');
       end
    %fluorescence should never be negative but sometimes with lots of 0s
    %looks like the smoothing f(x) makes it so
       SmoothParticles(n).Smoothed(SmoothParticles(n).Smoothed<0)=0;
       if exist('ShortTrace') & ShortTrace=='y'
       SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:end))]');
       else
       SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:(nc14+100)))]');
       end
%Slope thresholds
ONThreshold=377; 
OFFThreshold=-377; 
       
StartFramePool=[];
StartSlopePool=[];
EndSlopePool=[];
Starts=[find(SmoothParticles(n).RawSlope>=ONThreshold)];
Starts=min(Starts);
StartFramePool=[StartFramePool,Starts];
EndPts=[find(SmoothParticles(n).RawSlope<=OFFThreshold)];
EndPts=min(EndPts);
EndFramePool=[1000,EndPts];
if ~isempty(EndPts) & (~isempty(Starts))
while EndFramePool(end) ~= EndFramePool(end-1)
    for ss=EndFramePool(end):length(SmoothParticles(n).RawSlope)
        if SmoothParticles(n).RawSlope(ss) >=ONThreshold
            %if StartFramePool(end)~=ss
            StartFramePool=[StartFramePool,ss];
            StartSlopePool=[StartSlopePool, SmoothParticles(n).RawSlope(ss)];
            break
            %end
        end
    end
    for ss=StartFramePool(end):(length(SmoothParticles(n).RawSlope))
        if SmoothParticles(n).RawSlope(ss) <=OFFThreshold
            EndFramePool=[EndFramePool, ss];
            EndSlopePool=[EndSlopePool,SmoothParticles(n).RawSlope(ss)];
            break
        end
    end
    if sum(SmoothParticles(n).RawSlope(EndFramePool(end):end)>=ONThreshold)==0
%         EndFramePool(end+1)=EndFramePool(end);
%         StartFramePool(end+1)=StartFramePool(end);
        PotentialPeaks(n).Starts=[StartFramePool];
        PotentialPeaks(n).Ends=[EndFramePool(2:end)];
        break
    elseif sum(SmoothParticles(n).RawSlope(StartFramePool(end):end)<=OFFThreshold)==0
%         EndFramePool(end+1)=EndFramePool(end);
%         StartFramePool(end+1)=StartFramePool(end);
        PotentialPeaks(n).Starts=[StartFramePool];
        
        PotentialPeaks(n).Ends=[EndFramePool(2:end)];
        break
    end
end
if length(PotentialPeaks(n).Ends)>1 & (PotentialPeaks(n).Ends(end)==PotentialPeaks(n).Ends(end-1))
    PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(1:end-1)];
end
if length(PotentialPeaks(n).Starts)>1 & (PotentialPeaks(n).Starts(end)==PotentialPeaks(n).Starts(end-1))
    PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(1:end-1)];
end
BurstProperties(n).SmoothTrace=SmoothParticles(n).Smoothed;
%Get rid of repeated 1st and 2nd ON frames in case an OFF came before ON in
%original trace - added 7/27/18
if (length(PotentialPeaks(n).Starts)>=2)& (PotentialPeaks(n).Starts(1)==PotentialPeaks(n).Starts(2))
    PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(2:end)];
end
% If a trace starts with an OFF frame before an ON frame, need to ignore it
% 7/27/18
if PotentialPeaks(n).Ends(1) < PotentialPeaks(n).Starts(1)
    if length(PotentialPeaks(n).Ends)==1
        PotentialPeaks(n).Ends=nan;
        PotentialPeaks(n).Starts=nan;
    else
    PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(2:end)];
    end
end
if (isnan(PotentialPeaks(n).Ends)) & (isnan(PotentialPeaks(n).Starts))
    continue
end
%If start frames > end frames, ignore the last start since we can't say how
%long the burst lasts
if length(PotentialPeaks(n).Starts) > (length(PotentialPeaks(n).Ends))
    BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts(1:end-1));
    BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
    BurstProperties(n).InterBurst=ElapsedTime(PotentialPeaks(n).Starts(2:end-1))-ElapsedTime(PotentialPeaks(n).Ends(1:end-1));
    BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts(1:end-1)];
else
BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts);
BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
BurstProperties(n).InterBurst=ElapsedTime(PotentialPeaks(n).Starts(2:end))-ElapsedTime(PotentialPeaks(n).Ends(1:end-1));
BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts];
end

BurstProperties(n).OFFFrames=[PotentialPeaks(n).Ends];

if PotentialPeaks(n).Ends(1)==1
    BurstProperties(n).BurstAmplitude(1)=0;
else
    BurstProperties(n).BurstAmplitude(1)=max(SmoothParticles(n).Smoothed([1:PotentialPeaks(n).Ends(1)]));
end
    
    for pp=2:length(PotentialPeaks(n).Ends)
    BurstProperties(n).BurstAmplitude(pp)=max(SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(pp-1):PotentialPeaks(n).Ends(pp)]))%[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
    BurstProperties(n).BurstAmplitude(pp)=(BurstProperties(n).BurstAmplitude(pp))-(SmoothParticles(n).Smoothed(PotentialPeaks(n).Starts(pp)));
    end

BurstProperties(n).BurstAmplitude=BurstProperties(n).BurstAmplitude(BurstProperties(n).BurstAmplitude~=0); %get rid of 1st 0 bc above started at pp=2

%Burst size by integrating area under curve from ON point to ON point
for pp=1:(length(BurstProperties(n).ONFrames)-1)
    BurstProperties(n).BurstSize(pp)=(trapz([BurstProperties(n).ONFrames(pp):BurstProperties(n).ONFrames(pp+1)],BurstProperties(n).SmoothTrace(BurstProperties(n).ONFrames(pp):BurstProperties(n).ONFrames(pp+1))));
end
if isfield(BurstProperties, 'BurstSize')
BurstProperties(n).BurstSize(end+1)=(trapz([BurstProperties(n).ONFrames(end):BurstProperties(n).OFFFrames(end)], BurstProperties(n).SmoothTrace(BurstProperties(n).ONFrames(end):BurstProperties(n).OFFFrames(end))));
end
% if PotentialPeaks(n).Ends(1)==1
%     BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(2:end)]-1)]';
% else
% BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
% end
temp=find([CompiledParticles.Nucleus]==SpotDiff(n).Nucleus);
BurstProperties(n).TotalmRNAError=CompiledParticles(temp(1)).TotalmRNAError;

BurstProperties(n).FirstTimeOn=ElapsedTime(PotentialPeaks(n).Starts(1));
BurstProperties(n).NBursts=length(BurstProperties(n).BurstAmplitude);
%Frequency set btwn time of first ON and 50min into nc14 
if ((length(ElapsedTime))-nc14)<100
    BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(end)-ElapsedTime(PotentialPeaks(n).Starts(1)));
else
BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(nc14+100)-ElapsedTime(PotentialPeaks(n).Starts(1)));
end
BurstProperties(n).FrequencyAct=(BurstProperties(n).NBursts)/(ElapsedTime(PotentialPeaks(n).Ends(end))-ElapsedTime(PotentialPeaks(n).Starts(1)));
BurstProperties(n).APBin=SpotDiff(n).APBin;
BurstProperties(n).Nucleus=SpotDiff(n).Nucleus;
%BurstProperties(n).Interburst=
BurstProperties(n).TotalmRNA=SpotDiff(n).TotalmRNAOne;
BurstProperties(n).FractON=sum(BurstProperties(n).Duration)/(length(ElapsedTime(nc14:end))); 
BurstProperties(n).TotalElapsedTime=ElapsedTime(end)-ElapsedTime(nc14);

end
    else
        BurstProperties(n).APBin=SpotDiff(n).APBin;
        BurstProperties(n).TotalmRNA=nan;
    end
end


%% Do same for 2nd spots     
for n=1:length(SpotDiff)
    if (isfield(SpotDiff, 'SpotTwo')) & (length(SpotDiff(n).SpotTwo)>=3) & (length(SpotDiff(n).SpotTwo)==length(ElapsedTime))
       %try limiting frames to nc14 to reduce extreme smoothing
       if (length(ElapsedTime)-nc14)>=100
           SmoothParticles2(n).Smoothed=smooth(ElapsedTime(nc14:(nc14+100)),SpotDiff(n).SpotTwo(nc14:(nc14+100)),0.1,'lowess');
       else
           ShortTrace='y'
           SmoothParticles2(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiff(n).SpotTwo(nc14:end),0.1,'lowess');
       end
    %fluorescence should never be negative but sometimes with lots of 0s
    %looks like the smoothing f(x) makes it so
       SmoothParticles2(n).Smoothed(SmoothParticles2(n).Smoothed<0)=0;
       if exist('ShortTrace') & ShortTrace=='y';
          SmoothParticles2(n).RawSlope=diff(SmoothParticles2(n).Smoothed)./([diff(ElapsedTime(nc14:end))]');
       else
       SmoothParticles2(n).RawSlope=diff(SmoothParticles2(n).Smoothed)./([diff(ElapsedTime(nc14:(nc14+100)))]');
       end
%Slope thresholds
% ONThreshold=2000;
% OFFThreshold=-2000;
       %OFFThresholdSlope=dx<=Q;
% AboveLine=[SmoothParticles(ss).Smoothed];
% AboveLine(~ONThresholdSlope)=nan;
% AboveLine(OFFThresholdSlope)=nan;
% BelowLine=[SmoothParticles(ss).Smoothed];
% BelowLine(ONThresholdSlope)=nan;
% BelowLine(~OFFThresholdSlope)=nan;
StartFramePool=[];
StartSlopePool=[];
EndSlopePool=[];
Starts=[find(SmoothParticles2(n).RawSlope>=ONThreshold)];
Starts=min(Starts);
StartFramePool=[StartFramePool,Starts];
EndPts=[find(SmoothParticles2(n).RawSlope<=OFFThreshold)];
EndPts=min(EndPts);
EndFramePool=[1000,EndPts];
if ~isempty(EndPts) & (~isempty(Starts))
while EndFramePool(end) ~= EndFramePool(end-1)
    for ss=EndFramePool(end):length(SmoothParticles2(n).RawSlope)
        if SmoothParticles2(n).RawSlope(ss) >=ONThreshold
            %if StartFramePool(end)~=ss
            StartFramePool=[StartFramePool,ss];
            StartSlopePool=[StartSlopePool, SmoothParticles2(n).RawSlope(ss)];
            break
            %end
        end
    end
    for ss=StartFramePool(end):(length(SmoothParticles2(n).RawSlope))
        if SmoothParticles2(n).RawSlope(ss) <=OFFThreshold
            EndFramePool=[EndFramePool, ss];
            EndSlopePool=[EndSlopePool,SmoothParticles2(n).RawSlope(ss)];
            break
        end
    end
    if sum(SmoothParticles2(n).RawSlope(EndFramePool(end):end)>=ONThreshold)==0
%         EndFramePool(end+1)=EndFramePool(end);
%         StartFramePool(end+1)=StartFramePool(end);
        PotentialPeaks2(n).Starts=[StartFramePool];
        PotentialPeaks2(n).Ends=[EndFramePool(2:end)];
        break
    elseif sum(SmoothParticles2(n).RawSlope(StartFramePool(end):end)<=OFFThreshold)==0
%         EndFramePool(end+1)=EndFramePool(end);
%         StartFramePool(end+1)=StartFramePool(end);
        PotentialPeaks2(n).Starts=[StartFramePool];
        
        PotentialPeaks2(n).Ends=[EndFramePool(2:end)];
        break
    end
end
if length(PotentialPeaks2(n).Ends)>1 & (PotentialPeaks2(n).Ends(end)==PotentialPeaks2(n).Ends(end-1))
    PotentialPeaks2(n).Ends=[PotentialPeaks2(n).Ends(1:end-1)];
end
if length(PotentialPeaks2(n).Starts)>1 & (PotentialPeaks2(n).Starts(end)==PotentialPeaks2(n).Starts(end-1))
    PotentialPeaks2(n).Starts=[PotentialPeaks2(n).Starts(1:end-1)];
end
BurstProperties(n+length(SpotDiff)).SmoothTrace=SmoothParticles2(n).Smoothed;
%Get rid of repeated 1st and 2nd ON frames in case an OFF came before ON in
%original trace - added 7/27/18
if (length(PotentialPeaks2(n).Starts)>=2) & (PotentialPeaks2(n).Starts(1)==PotentialPeaks2(n).Starts(2))
    PotentialPeaks2(n).Starts=[PotentialPeaks2(n).Starts(2:end)];
end
% If a trace starts with an OFF frame before an ON frame, need to ignore it
% 7/27/18
if PotentialPeaks2(n).Ends(1) < PotentialPeaks2(n).Starts(1)
    if length(PotentialPeaks2(n).Ends)==1
        PotentialPeaks2(n).Ends=nan;
        PotentialPeaks2(n).Starts=nan;
    else
    PotentialPeaks2(n).Ends=[PotentialPeaks2(n).Ends(2:end)];
    end
end
%if just said ends and starts are nan bc only 1 of each that not same
%burst, jump to next spot
if (isnan(PotentialPeaks2(n).Ends)) & (isnan(PotentialPeaks2(n).Starts))
    continue
end
%If start frames > end frames, ignore the last start since we can't say how
%long the burst lasts
if length(PotentialPeaks2(n).Starts) > (length(PotentialPeaks2(n).Ends))
    BurstProperties(n+length(SpotDiff)).Duration=ElapsedTime(PotentialPeaks2(n).Ends)-ElapsedTime(PotentialPeaks2(n).Starts(1:end-1));
    BurstProperties(n+length(SpotDiff)).Duration=BurstProperties(n+length(SpotDiff)).Duration(BurstProperties(n+length(SpotDiff)).Duration>=0);
    BurstProperties(n+length(SpotDiff)).InterBurst=ElapsedTime(PotentialPeaks2(n).Starts(2:end-1))-ElapsedTime(PotentialPeaks2(n).Ends(1:end-1));
    BurstProperties(n+length(SpotDiff)).ONFrames=[PotentialPeaks2(n).Starts(1:end-1)];
else
BurstProperties(n+length(SpotDiff)).Duration=ElapsedTime(PotentialPeaks2(n).Ends)-ElapsedTime(PotentialPeaks2(n).Starts);
BurstProperties(n+length(SpotDiff)).Duration=BurstProperties(n+length(SpotDiff)).Duration(BurstProperties(n+length(SpotDiff)).Duration>=0);
BurstProperties(n+length(SpotDiff)).InterBurst=ElapsedTime(PotentialPeaks2(n).Starts(2:end))-ElapsedTime(PotentialPeaks2(n).Ends(1:end-1));
end
BurstProperties(n+length(SpotDiff)).ONFrames=[PotentialPeaks2(n).Starts];
BurstProperties(n+length(SpotDiff)).OFFFrames=[PotentialPeaks2(n).Ends];
if PotentialPeaks2(n).Ends(1)==1
    BurstProperties(n+length(SpotDiff)).BurstAmplitude(1)=0;
else
    BurstProperties(n+length(SpotDiff)).BurstAmplitude(1)=max(SmoothParticles2(n).Smoothed([1:PotentialPeaks2(n).Ends(1)]));
end
 %if length(PotentialPeaks2(n).Ends)>=2   
    for pp=2:length(PotentialPeaks2(n).Ends)
    BurstProperties(n+length(SpotDiff)).BurstAmplitude(pp)=max(SmoothParticles2(n).Smoothed([PotentialPeaks2(n).Ends(pp-1):PotentialPeaks2(n).Ends(pp)]))%[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
    BurstProperties(n+length(SpotDiff)).BurstAmplitude(pp)=(BurstProperties(n+length(SpotDiff)).BurstAmplitude(pp))-(SmoothParticles2(n).Smoothed(PotentialPeaks2(n).Starts(pp)));
    end
 %end

BurstProperties(n+length(SpotDiff)).BurstAmplitude=BurstProperties(n+length(SpotDiff)).BurstAmplitude(BurstProperties(n+length(SpotDiff)).BurstAmplitude~=0); %get rid of 1st 0 bc above started at pp=2

%Burst size by integrating area under curve from ON point to ON point
for pp=1:(length(BurstProperties(n+length(SpotDiff)).ONFrames)-1)
    BurstProperties(n+length(SpotDiff)).BurstSize(pp)=(trapz([BurstProperties(n+length(SpotDiff)).ONFrames(pp):BurstProperties(n+length(SpotDiff)).ONFrames(pp+1)],BurstProperties(n+length(SpotDiff)).SmoothTrace(BurstProperties(n+length(SpotDiff)).ONFrames(pp):BurstProperties(n+length(SpotDiff)).ONFrames(pp+1))));
end
BurstProperties(n+length(SpotDiff)).BurstSize(end+1)=(trapz([BurstProperties(n+length(SpotDiff)).ONFrames(end):BurstProperties(n+length(SpotDiff)).OFFFrames(end)], BurstProperties(n+length(SpotDiff)).SmoothTrace(BurstProperties(n+length(SpotDiff)).ONFrames(end):BurstProperties(n+length(SpotDiff)).OFFFrames(end))));

% if PotentialPeaks2(n).Ends(1)==1
%     BurstProperties(n+length(SpotDiff)).BurstAmplitude=[SmoothParticles2(n).Smoothed([PotentialPeaks2(n).Ends(2:end)]-1)]';
% else
% BurstProperties(n+length(SpotDiff)).BurstAmplitude=[SmoothParticles2(n).Smoothed([PotentialPeaks2(n).Ends]-1)]';
% end
temp2=find([CompiledParticles.Nucleus]==SpotDiff(n).Nucleus);
BurstProperties(n+length(SpotDiff)).TotalmRNAError=CompiledParticles(temp2(2)).TotalmRNAError;

BurstProperties(n+length(SpotDiff)).FirstTimeOn=ElapsedTime(PotentialPeaks2(n).Starts(1));
BurstProperties(n+length(SpotDiff)).NBursts=length(BurstProperties(n+length(SpotDiff)).BurstAmplitude);
% frequency set btwn time of first burst and 50min into nc14
if ((length(ElapsedTime))-nc14)<100
    BurstProperties(n+length(SpotDiff)).Frequency=(BurstProperties(n+length(SpotDiff)).NBursts)/(ElapsedTime(end)-ElapsedTime(PotentialPeaks2(n).Starts(1)));
else
BurstProperties(n+length(SpotDiff)).Frequency=(BurstProperties(n+length(SpotDiff)).NBursts)/(ElapsedTime(nc14+100)-ElapsedTime(PotentialPeaks2(n).Starts(1)));
end
BurstProperties(n+length(SpotDiff)).FrequencyAct=(BurstProperties(n+length(SpotDiff)).NBursts)/(ElapsedTime(PotentialPeaks2(n).Ends(end))-ElapsedTime(PotentialPeaks2(n).Starts(1)));
BurstProperties(n+length(SpotDiff)).APBin=SpotDiff(n).APBin;
BurstProperties(n+length(SpotDiff)).Nucleus=SpotDiff(n).Nucleus;
BurstProperties(n+length(SpotDiff)).TotalmRNA=SpotDiff(n).TotalmRNATwo;
%BurstProperties(n+length(SpotDiff)).TotalmRNAError=SpotDiff(n).
BurstProperties(n+length(SpotDiff)).FractON=sum(BurstProperties(n+length(SpotDiff)).Duration)/length(ElapsedTime(nc14:end));
BurstProperties(n+length(SpotDiff)).TotalElapsedTime=ElapsedTime(end)-ElapsedTime(nc14);
end
    else
        BurstProperties(n+length(SpotDiff)).APBin=SpotDiff(n).APBin;
        BurstProperties(n+length(SpotDiff)).TotalmRNA=nan;
    end
    
end

for bb=1:length(BurstProperties)
    if isempty(BurstProperties(bb).APBin)
        BurstProperties(bb).APBin=nan;
        BurstProperties(bb).Nucleus=nan;
    end
   
end
clear ShortTrace
%save('BurstPropertiesSlope','BurstProperties');