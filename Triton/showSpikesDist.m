function [inDist,inSort,resDist,resSort]=showSpikesDist(Obj,resSpikeList,resSpikeLegend,show) % Input is DATA structure
%resSpikeList
if(nargin<1); resSpikeList ={'RES'};end
if(nargin<4); show = 0 ;end;
if(nargout<1); show = 1; end;
if(nargin<3 || numel(resSpikeLegend)<numel(resSpikeList)); resSpikeLegend=resSpikeList; end

[inputSpikes,reservoirSpikes,inDist,resDist] = getDATAres(Obj,resSpikeList);
[inDist,inSort] = sort(inDist);
[resDist,resSort] =sort(resDist);
if(show)
    hold on;
    h0=fill([1 1:length(inDist) length(inDist)],[0; inDist; 0],'r','FaceAlpha',0.2,'EdgeAlpha',0.1);
    line([  0 length(inDist)]',repmat([mean(inDist)],1,2)','Color','r','LineStyle','--','LineWidth',2);
    distLeg{1} = ['Input (#:' num2str(inputSpikes,'%.0f') ' s^{-1})'];
    Color = ['g','b','y','m'];
    h=[];
    for iRes = 1:numel(resSpikeList)
        h(iRes)=fill([1 1:length(resDist(:,iRes)) length(resDist(:,iRes))],[0; resDist(:,iRes); 0],Color(mod(iRes-1,length(Color))+1),'FaceAlpha',0.2,'EdgeAlpha',0.1);
        distLeg{iRes+1} = [resSpikeLegend{iRes} '(#:' num2str(reservoirSpikes(iRes),'%.0f') ' s^{-1})'];
        line([ 0 length(resDist(:,iRes))]',repmat([mean(resDist(:,iRes))],1,2)','Color',Color(mod(iRes-1,length(Color))+1),'LineStyle','--','LineWidth',2);
    end
    %line([190 190;125 125]',[0 max(inputDist);0 max(reservoirDist)]','LineStyle','--');
    
    legend([h0 h],distLeg);
    xlabel('Neuron(after sorting)'); ylabel('spike rate(s^{-1}/neuron)'); title('Spikes Rate Distribution');
end
end
function [iS,rS,inS,rnS]=getDATAres(Obj,resSpikeList)
[~,simTime] = size([Obj(:).S]);
epochTime = simTime*1E-3;
% return #inputSpikes (iS) and # reservoir Spikes (rS), n for normalized
inS = sum([Obj(:).S],2)/epochTime;
iS = sum(inS,1);
rnS = [];
for iRes = 1:numel(resSpikeList)
    rnS(:,iRes) = sum([Obj(:).(resSpikeList{iRes})],2)/epochTime;
end
rS = sum(rnS,1);
end
function compareDist(Obj)

[iS,rS,inS,rnS]=getDATAres(Obj);hold on; plot(rS);
end
