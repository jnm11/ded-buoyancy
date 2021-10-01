function typ=ded_gc_f7_g_classification(nm)
% typ=ded_gc_f7_g_classification(nm) return instability type
% S stable
% H Holmboe
% K Kelvin Helmholtz

if iscell(nm)
  for j=1:length(nm)
    typ{j}=ded_gc_f7_g_classification(nm{j});
  end
  return;
end

typ=[];
a=cellstrtok(nm,'/');
for j=1:length(a)
  if strcmp(a{j},'gc')
    break;
  end
end
a={a{j:end}};
if isempty(a) 
  error(sprintf('ded_gc_f7_g_classification: %s is not a gravity current simulation',nm));
end

switch(a{2})
  case({'emle','ccle'})
    typ='K';
  case('le')
    switch(a{3})
      case('a')
        Pe=str2num(a{4});        Re=str2num(a{5});        Nx=str2num(a{6});
        switch(Pe)
          case(4000);typ=shk(Re,[250 400    ], {'S','M','K'});
          case(8000);typ=shk(Re,[250 400 600 ],{'S','H','M','K'}); % 'M' mixed Holmboe KH
        end
      case('a')
        Pe=str2num(a{4});        Re=str2num(a{5});        Nx=str2num(a{6});
        switch(Pe)
          case(4000);typ=shk(Re,[400 500    ], {'S','M','K'});
          case(8000);typ=shk(Re,[250 400 600 ],{'S','H','M','K'}); % 'M' mixed Holmboe KH
        end
      end  
  case('f6')
    switch(a{3})
      case('f')
        typ='K';
      case({'g','h'})
        Pe=str2num(a{4}); Re=str2num(a{5}); Nx=str2num(a{6}); g=str2num(a{7});
        switch(Pe)
          case(1000); typ=shk(Re,[1200     ],{'S','K'});
          case(2000); typ=shk(Re,[1200     ],{'S','K'});
          case(4000); typ=shk(Re,[0700     ],{'S','K'});
          case(8000); typ=shk(Re,[0500 0800],{'S','H','K'});
        end
      case({'i','mg'})
        Pe=str2num(a{4}); Re=str2num(a{5}); Nx=str2num(a{6}); g=str2num(a{7});
        switch(Pe)
          case(1000); typ=shk(Re,[1200     ],{'S','K'});
          case(2000); typ=shk(Re,[1200     ],{'S','K'});
          case(4000); typ=shk(Re,[0700     ],{'S','K'});
          case(8000); typ=shk(Re,[0400 0800],{'S','H','K'});
        end
    end
end

if isempty(typ)
  keyboard;
end


function typ=shk(Re,RR,tt)
typ=tt{end};
for k=1:length(RR)
  if Re<RR(k)
    typ=tt{k};
    break;
  end
end

