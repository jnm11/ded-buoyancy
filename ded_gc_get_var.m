function [x,xs]=ded_gc_get_var(c,nm)

xs=NaN;

switch(nm)
  case('Nx'); x=c.Nx;
  case('Ny'); x=c.Ny;
  case('Nz'); x=c.Nz;
  case('g');  x=c.g;if isfield(c,'sg');xs=c.sg;end
  case('X');  x=c.X;if isfield(c,'sX');xs=c.sX;end
  case('U');  x=c.U;if isfield(c,'sU');xs=C.sU;end
  case('L');  x=c.L;
  case('H');  x=c.H;
  case('W');  x=c.W;
  case('dx'); x=round(c.L./c.Nx*1e6)/1e6;
  case('R');  x=c.R;
  case('Re'); x=c.Re;
  case('Scb');x=c.Scb;
  case('Scs');x=c.Scs;
  case('Peb');x=c.Scb.*c.Re;
  case('Pes');x=c.Scs.*c.Re;
  case('sg'); x=c.sg;
  case('sX'); x=c.sX;
  case('sR'); x=c.sR;
  case('H');  x=c.H;
  case('hu'); x=c.hu;
  case('hb'); x=c.hb;
  case('U');  x=c.U;
  case('nm1');x=c.nm1;
  case('nm2');x=c.nm2;
  case('nm3');x=c.nm3;
  case('nm4');x=c.nm4;
  case('nm5');x=c.nm5;
  case('nma');x=c.nma;
  case('nmb');x=c.nmb;
  case('nmc');x=c.nmc;
  case('nmd');x=c.nmd;
  case('nme');x=c.nme;
end
