function a=ded_read_flux(n)

a1=ded_read_g(n,'flux');
a2=ded_read_g(n,'yz');
if isempty(a1) & isempty(a2)
  a=[];
  return;
elseif ~isempty(a1) & ~isempty(a2)
  a1=lower_fieldnames(a1);
  a.t  = [    a1.t; a2.t]; 
  a.b  = [    a1.b  a2.b]; 
  a.p  = [    a1.p  a2.p]; 
  a.bu = [    a1.q  a2.bu]; 
  a.uu = [   a1.uu  a2.uu];
  a.uw = [   a1.uw  a2.uw]; 
  a.vv = [   a1.vv  a2.vv]; 
  a.ww = [   a1.ww  a2.ww]; 
  
  a.x=a2.x;

  a.E  = [ NaN*a1.uu a2.E  ];
  a.S  = [ NaN*a1.uu a2.S  ];
  a.bb = [ NaN*a1.uu a2.bb ];
  a.bu = [ NaN*a1.uu a2.bu ];
  a.bv = [ NaN*a1.uu a2.bv ];
  a.bw = [ NaN*a1.uu a2.bw ];
  a.u  = [ NaN*a1.uu a2.u  ];
  a.uv = [ NaN*a1.uu a2.uv ];
  a.v  = [ NaN*a1.uu a2.v  ];
  a.vw = [ NaN*a1.uu a2.vw ];
  a.w  = [ NaN*a1.uu a2.w  ];
elseif ~isempty(a1) 
  a=lower_fieldnames(a1);
else
  a=a2;
end
if isempty(a)
  return;
end
a=lower_fieldnames(a);
f=findstr('dedalus',n);
if isempty(f)
  a.nm=n;
else
  a.nm=n(f+8:end);
end


  
