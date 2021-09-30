function ded_print(n)
if nargin<1
  n=[];
end
if isempty(n)
  n='gc/*';
end

[a A]=ded_all_stats(n);
if isempty(a)
  return;
end
ded_stats_plot(a);
B=A;
FR = A(:,1); % Forcing
Nx = A(:,2); % x grid size
L  = A(:,3); % simulation length
H  = A(:,4); % simulation height
W  = A(:,5); % simulation width
U  = A(:,6);
U1 = A(:,7);
Re = A(:,8); % simulation width
PIDX = A(:,9); % simulation width

A=A(:,1:5);
A0=0*A(1,:);

fp=fopen('~/misc/gc-status','w');
n=length(a);
for j=1:n
  if any(A(j,:)~=A0)
    e=sprintf('\n\nf=%1i, Nx=%4i, Re=%5.0f, PX=%5.0f, L=%5.1f, H=%3.1f, W=%3.1f\n',FR(j),Nx(j),Re(j),PIDX(j),L(j),H(j),W(j));
    fprintf(1,e);
    fprintf(fp,e);
  end
  A0=A(j,:);
  e=sprintf('f=%1i, Nx=%4i, Re=%5.0f, PX=%5.0f, L=%5.1f, H=%3.1f, W=%3.1f, U=%6.3f, %s %s\n',FR(j),Nx(j),Re(j),PIDX(j),L(j),H(j),W(j),U(j),a(j).msg,a(j).status);
  fprintf(1,e);
  fprintf(fp,e);
end
fclose(fp);
unix('cat ~/misc/gc-status');

return;
