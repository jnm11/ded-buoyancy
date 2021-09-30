function [b a]=ded_read_state(fn,p)
%a=ded_read_state(fn) Read a final or checkpoint file
%fn='gc/qgc2/01/checkpoint/checkpoint_s2.hdf5';

if isfile(fn)
  a=ded_read_hdf(fn);
else
  fn1=[ded_dedalus_data_dir '/' fn];
  if isfile(fn)
    a=ded_read_hdf(fn1);
  else
    disp(sprintf('ded_read_state: could not find a file %s',fn));
    a=[];
    b=[];
    return
  end
end

nms=fieldnames(a);
if isfield(a,'Tz')
  chb=1;
else
  chb=0;
end

if p.Ny>1
  n=[p.Nz p.Ny p.Nx];
else
  n=[p.Nz p.Nx];
end


for j=1:length(nms)
  nm=nms{j};
  if isstruct(a.(nm))
    z=complex(a.(nm).r,a.(nm).i);
    a.(nm)=z;
    if ndims(z)>2
      if n(3)==2*size(z,3) 
        z=ifft(cat(3,z,zeros([size(z,1) size(z,2) 1]),conj(z(:,:,end:-1:2))),[],3)*n(3);
      else
        keyboard
      end
    end 
    if ndims(z)>1
      if n(2)==2*size(z,2) 
        z=ifft(cat(2,z,zeros([size(z,1) 1 size(z,3)]),conj(z(:,end:-1:2,:))),[],2)*n(2);
      elseif n(2)==size(z,2)+1
        z=ifft(cat(2,z(:,1:n(2)/2,:),zeros([size(z,1) 1 size(z,3)]),conj(z(:,n(2)/2:-1:2,:))),[],2)*n(2);
      elseif n(2)==size(z,2) % sincos basis should actually take account of parity
        z=ifft(cat(2,z(:,1:n(2)/2,:),conj(z(:,n(2)/2:-1:2,:))),[],2)*n(2);
      else
        keyboard
      end
    end
    if chb
      z=ichebc2f(z,1,[],[0 p.H]); % This messes up complex numbers
    else
      if n(1)==2*size(z,1) 
        z=ifft(cat(1,z,zeros([1 size(z,2) size(z,3)]),conj(z(end:-1:2,:,:))),[],1)*n(1);
      elseif n(1)==size(z,1)+1
        z=ifft(cat(1,z(1:n(1)/2,:,:),zeros([1 size(z,2) size(z,3)]),conj(z(n(1)/2:-1:2,:,:))),[],1)*n(1);
      else
        keyboard
      end
    end
    e=max(abs(imag(z(:))));
    if e>1e-5
      disp(sprintf('ded_read_state: Field %s is not real %8.2e',nm,e));
    end
    b.(nm)=real(z);              
  else
    b.(nm)=a.(nm);
  end
end
if chb
  b.z=p.H/2*(1+ichebgrid(size(z,1)));
end
return
 
max(abs(imag(z(:))))


nm=[ded_dedalus_data_dir '/gc/slice/00'];
sfn=[nm '/checkpoint/checkpoint_s12.hdf5'];
sfn=[nm '/final/final_s568.hdf5'];
ufn=[nm '/u/u_s568.hdf5'];
bfn=[nm '/b/b_s568.hdf5'];
ffn=[nm '/force/force_s1.hdf5'];
p=ded_read_param(nm);
[a c]=ded_read_state(sfn,p);
u=ded_read_hdf(ufn);
b=ded_read_hdf(bfn);
f=ded_read_hdf(ffn);
%f=ded_zgrid(f,[],[],{},{},{},p.H);
%u=ded_zgrid(u,[],[],{},{},{},p.H);
%c=ded_zgrid(u,[],[],{},{},{},p.H);
a.x=u.x;
c.x=u.x;
%a=ded_zgrid(a,2*length(a.Tz),[],{},{},p.H);

mesh(a.x,a.z,a.u-u.u);view([50 60]);


subplot(2,2,1);
mesh(u.x,u.z,u.u);view([50 60]);
subplot(2,2,3);
mesh(f.x,f.z,f.wu);view([50 60]);
subplot(2,2,2);
mesh(f.x,f.z,f.fu);view([50 60]);
subplot(2,2,4);
mesh(a.x,a.z,a.u);view([50 60]);

A=ichebc2f(c.u,1,[],[0 p.H]);
a.u=ifft([A 0*A(:,1) conj(A(:,end:-1:2))],[],2,'symmetric');

u=ifft(A,[],2,'symmetric');



z=[a.u 0*a.u(:,1) conj(a.u(:,end:-1:2))];
mesh(ifft(z,[],2,'symmetric'));
z=ichebc2f(z,1);
mesh(u.x,u.z,ifftc.u);view([50 60]);

f=ded_read_hdf(ffn);
f=ded_zgrid(f,[],[],[],[],1,p.H);


sum(f.wu(:))*(f.x(2)-f.x(1))*(f.z(2)-f.z(1))/(p.L*p.H*p.W)

wuI=sum(f.wuI(end,:))*(f.x(2)-f.x(1))/(p.L*p.H*p.W);
wbI=sum(f.wbI(end,:))*(f.x(2)-f.x(1))/(p.L*p.H*p.W);
disp(sprintf('Mean wu:%8.4f wb:%8.4f',wuI,wbI'));
0.041021675442836275
0.08789283905368495
