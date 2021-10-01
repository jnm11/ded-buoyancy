function p=ded_force_fx(nm,trg)
if nargin<2
  trg=[0 inf];
end
c=ded_coord(nm);
p=ded_read_param(nm);
w=ded_hdf(nm,'fx','force');
if isempty(w); return ; end;

if isfield(w,'sdl');nsdl=length(w.sdl);x.sdl=c.Ax(1:nsdl);Isdl=sum(w.sdl*c.dAx);s.sdl=sprintf('left d sample %4.2f',Isdl);end
if isfield(w,'wdl');nwdl=length(w.wdl);x.wdl=c.Ax(1:nwdl);Iwdl=sum(w.wdl*c.dAx);s.wdl=sprintf('left d source %4.2f',Iwdl);end
if isfield(w,'sul');nsul=length(w.wul);x.sul=c.Ax(1:nsul);Isul=sum(w.sul*c.dAx);s.sul=sprintf('left u sample %4.2f',Iwul);end
if isfield(w,'wul');nwul=length(w.wul);x.wul=c.Ax(1:nwul);Iwul=sum(w.wul*c.dAx);s.wul=sprintf('left u source %4.2f',Iwul);end
if isfield(w,'sbl');nsbl=length(w.sbl);x.sbl=c.Ax(1:nsbl);Isbl=sum(w.sbl*c.dAx);s.sbl=sprintf('left b sample %4.2f',Isbl);end
if isfield(w,'wbl');nwbl=length(w.wbl);x.wbl=c.Ax(1:nwbl);Iwbl=sum(w.wbl*c.dAx);s.wbl=sprintf('left b source %4.2f',Iwbl);end

if isfield(w,'sdr');nsdr=length(w.sdr);x.sdr=c.Ax(end-nsdr+1:end);Isdr=sum(w.sdr*c.dAx);s.sdr=sprintf('left d sample %4.2f',Isdr);end
if isfield(w,'wdr');nwdr=length(w.wdr);x.wdr=c.Ax(end-nwdr+1:end);Iwdr=sum(w.wdr*c.dAx);s.wdr=sprintf('left d source %4.2f',Iwdr);end
if isfield(w,'sur');nsur=length(w.wur);x.sur=c.Ax(end-nsur+1:end);Isur=sum(w.sur*c.dAx);s.sur=sprintf('left u sample %4.2f',Iwur);end
if isfield(w,'wur');nwur=length(w.wur);x.wur=c.Ax(end-nwur+1:end);Iwur=sum(w.wur*c.dAx);s.wur=sprintf('left u source %4.2f',Iwur);end
if isfield(w,'sbr');nsbr=length(w.sbr);x.sbr=c.Ax(end-nsbr+1:end);Isbr=sum(w.sbr*c.dAx);s.sbr=sprintf('left b sample %4.2f',Isbr);end
if isfield(w,'wbr');nwbr=length(w.wbr);x.wbr=c.Ax(end-nwbr+1:end);Iwbr=sum(w.wbr*c.dAx);s.wbr=sprintf('left b source %4.2f',Iwbr);end

lrnm{1}=intersect(fieldnames(w),{'sdl' ,'wdl' ,'sul' ,'wul' ,'sbl' ,'wbl'});                 
lrnm{2}=intersect(fieldnames(w),{'sdr' ,'wdr' ,'sur' ,'wur' ,'sbr' ,'wbr'});                 
xs={'x0','x1','x2','x3','x4','x5','x6','x7'};
figure(1);clf;
co=get(0,'defaultaxescolororder');

for k=1:2
  subplot(2,1,k);
  hold('on');h=[];
  for j=1:length(lrnm{k})
    fnm=lrnm{k}{j};
    h(j)=plot(x.(fnm),w.(fnm)/max(w.(fnm)),'-s','color',co(j,:));
  end
  for j=1:length(xs)
    if isfield(p,xs{j})
      line([1;1]*p.xs{j},'color',0.7*[1 1 1]);
    end
  end
  legend(h,lrnm{k},'location','best');
  axis([-inf inf -inf inf])
end


