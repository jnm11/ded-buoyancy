function nms=ded_dir_nms(nm) 
DD=[ded_dedalus_data_dir '/'];
nms=cellstr_ls([DD nm],'ls -d');
