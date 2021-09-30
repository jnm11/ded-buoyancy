function p=ded_read_parity(nm)
fnp=[ded_dedalus_data_dir '/' nm '/' 'parity.hdf5'];
if isfile(fnp)
  p=ded_read_hdf([ded_dedalus_data_dir '/' nm '/' 'parity.hdf5']);
  return;
end
fnp=[ded_dedalus_data_dir '/results/' nm '/' 'parity.mat'];
if isfile(fnp)
  p=load(fnp);
  p=p.p;
  return;
end

a = ded_read_param(nm);
switch(a.sType)
  case 'gc'
    if strcmp(a.Tx,'SinCos'); 
      ux = -1;
      vx =  1;
      wx =  1;
      px =  1;
      bx =  1;
    else
      ux = 0;
      vx = 0;
      wx = 0;
      px = 0;
      bx = 0;
    end
    if isfield(a,'Ty')
      if strcmp(a.Ty,'SinCos'); 
        uy =  1;
        vy = -1;
        wy =  1;
        py =  1;
        by =  1;
      else
        uy =  0;
        vy =  0;
        wy =  0;
        py =  0;
        by =  0;
      end
    else
      uy =  [];
      vy =  [];
      wy =  [];
      py =  [];
      by =  [];
    end
    if strcmp(a.Tz,'SinCos'); 
      uz =  1;
      vz = -1;
      wz =  1;
      pz =  1;
      bz =  1;
    else
      uz =  0;
      vz =  0;
      wz =  0;
      pz =  0;
      bz =  0;
    end
end
p.u=[ux uy uz];
p.v=[vx vy vz];
p.w=[wx wy wz];
p.p=[px py pz];
p.b=[bx by bz];

nm={'u','v','w','p','b'};
for j=length(nm)
  for k=1:length(nm)
    p.([nm{j} nm{k}])=p.(nm{j}).*p.(nm{k});
    end
end
nm=fieldnames(p);

if a.Ny==1
  for j=length(nm)
    p.(['d' nm{j} 'dx'])= p.(nm{j}).*[-1 +1];
    p.(['d' nm{j} 'dz'])= p.(nm{j}).*[+1 -1];
  end
  p.c='xz';
  nm=fieldnames(p);
  for j=length(nm)
    if any(nm{j}=='y') | any(nm{j}=='v')
      p=rmfield(p,nm{j});
    else
      p.(nm{j})=nm{j}([1 3]);
    end
  end
else
  for j=length(nm)
    p.(['d' nm{j} 'dx'])= p.(nm{j}).*[-1 +1 +1];
    p.(['d' nm{j} 'dy'])= p.(nm{j}).*[+1 -1 +1];
    p.(['d' nm{j} 'dz'])= p.(nm{j}).*[+1 +1 -1];
  end
  p.c='xyz';
end


      
