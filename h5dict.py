import h5py
import numpy as np
from pathlib import Path

# Write parameters to a new file then rename the old one, then rename the new one
def save_dict_to_hdf5(dic, filename):
    sfn=str(filename)
    fnew=sfn + '.new'
    
    with h5py.File(fnew, 'w') as h5file:
        fn=Path(sfn)
        fnew=Path(fnew)
        recursively_save_dict_contents_to_group(h5file, '/', dic)
        if fn.is_file():   fn.rename(sfn + '.old')
        if fnew.is_file(): fnew.rename(fn)
        
def load_dict_from_hdf5(filename):
    with h5py.File(str(filename), 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')



def recursively_save_dict_contents_to_group( h5file, path, dic):
    
    # argument type checking
    if not isinstance(dic, dict):
        raise ValueError("must provide a dictionary")        

    if not isinstance(path, str):
        raise ValueError("path must be a string")
    if not isinstance(h5file, h5py._hl.files.File):
        raise ValueError("must be an open h5py file")
    # save items to the hdf5 file
    for key, item in dic.items():
        #print(key,item)
        key = str(key)
        if isinstance(item, list):
            item = np.array(item)
            #print(item)
        if not isinstance(key, str):
            raise ValueError("dict keys must be strings to save to hdf5")
        # save strings, numpy.int64, and numpy.float64 types
        if isinstance(item, (np.int64, np.float64, str, np.float, float, np.float32,int)):
            #print( 'here' )
            h5file[path + key] = item
#            if not h5file[path + key][()] == item:
#                raise ValueError('The data representation in the HDF5 file does not match the original dict.')
        elif isinstance(item, (np.bool,np.bool_)):
            #print( 'here' )
            if item: item=np.int(1)
            else:    item=np.int(0)
            h5file[path + key] = item
#            if not h5file[path + key][()] == item:
#                raise ValueError('The data representation in the HDF5 file does not match the original dict.')
        # save numpy arrays
        elif isinstance(item, np.ndarray):            
            try:
                h5file[path + key] = item
            except:
                item = np.array(item).astype('|S9')
                h5file[path + key] = item
#            if not np.array_equal(h5file[path + key][()], item):
#                raise ValueError('The data representation in the HDF5 file does not match the original dict.')
        # save dictionaries
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        # other types cannot be saved and will result in an error
        else:
            #print(item)
            raise ValueError('Cannot save %s type.' % type(item))

def recursively_load_dict_contents_from_group( h5file, path): 

    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans         
