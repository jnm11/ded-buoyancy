function d=ded_dedalus_data_dir()
d=getenv('DEDALUS_DATA');
if isempty(d)
  error('ded_read_param: then environment variable DEDALUS_DATA must be set');
end
%$DEDALUS_DATA
