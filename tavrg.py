
def init():
    self.sz=[Nx,Ny]

from mpi4py import MPI
import h5py


class tavrg(Handler):
    """
    Handler that time averages fields in a memory efficient manner
    taken from FileHandler class in evaluator

    rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)
    f = h5py.File('parallel_test.hdf5', 'w', driver='mpio', comm=MPI.COMM_WORLD)
    dset = f.create_dataset('test', (4,), dtype='i')
    dset[rank] = rank
    f.close()
    

     Parameters
    ----------
    base_path : str
        Base path for analyis output folder
    max_writes : int, optional
        Maximum number of writes per set (default: infinite)
    max_size : int, optional
        Maximum file size to write to, in bytes (default: 2**30 = 1 GB).
        (Note: files may be larger after final write.)
    parallel : bool, optional
        Perform parallel writes from each process to single file (True), or
        separately write to individual process files (False).
        Default behavior set by config option.
    mode : str, optional
        'overwrite' to delete any present analysis output with the same base path.
        'append' to begin with set number incremented past any present analysis output.
        Default behavior set by config option.

    """
     def __init__(self, base_path, *args, max_writes=np.inf, max_size=2**30, mode=None, **kw):

        Handler.__init__(self, *args, **kw)

        if mode is None: mode = FILEHANDLER_MODE_DEFAULT

        # Check base_path
        base_path = pathlib.Path(base_path).absolute()
        if any(base_path.suffixes):
            raise ValueError("base_path should indicate a folder for storing HDF5 files.")

        # Attributes
        self.base_path = base_path
        self.max_writes = max_writes
        self.max_size = max_size
        self.parallel = parallel
        self._sl_array = np.zeros(1, dtype=int)

        # Resolve mode
        mode = mode.lower()
        if mode not in ['overwrite', 'append']:
            raise ValueError("Write mode {} not defined.".format(mode))

        comm = self.domain.dist.comm_cart
        if comm.rank == 0:
            set_pattern = '%s_s*' % (self.base_path.stem)
            sets = list(self.base_path.glob(set_pattern))
            if mode == "overwrite":
                for set in sets:
                    if set.is_dir():
                        shutil.rmtree(str(set))
                    else:
                        set.unlink()
                set_num = 1
                total_write_num = 1
            elif mode == "append":
                set_nums = []
                if sets:
                    for set in sets:
                        m = re.match("{}_s(\d+)$".format(base_path.stem), set.stem)
                        if m:
                            set_nums.append(int(m.groups()[0]))
                    max_set = max(set_nums)
                    joined_file = base_path.joinpath("{}_s{}.h5".format(base_path.stem,max_set))
                    p0_file = base_path.joinpath("{0}_s{1}/{0}_s{1}_p0.h5".format(base_path.stem,max_set))
                    if os.path.exists(str(joined_file)):
                        with h5py.File(str(joined_file),'r') as testfile:
                            last_write_num = testfile['/scales/write_number'][-1]
                    elif os.path.exists(str(p0_file)):
                        with h5py.File(str(p0_file),'r') as testfile:
                            last_write_num = testfile['/scales/write_number'][-1]
                    else:
                        last_write_num = 0
                        logger.warn("Cannot determine write num from files. Restarting count.")
                else:
                    max_set = 0
                    last_write_num = 0
                set_num = max_set + 1
                total_write_num = last_write_num + 1
        else:
            set_num = None
            total_write_num = None
        # Communicate set and write numbers
        self.set_num = comm.bcast(set_num, root=0)
        self.total_write_num = comm.bcast(total_write_num, root=0)

        # Create output folder
        with Sync(comm):
            if comm.rank == 0:
                base_path.mkdir(exist_ok=True)

        # Set HDF5 property list for collective writing
        self._property_list = h5py.h5p.create(h5py.h5p.DATASET_XFER)
        self._property_list.set_dxpl_mpio(h5py.h5fd.MPIO_COLLECTIVE)

    def check_file_limits(self):
        """Check if write or size limits have been reached."""
        write_limit = (self.file_write_num >= self.max_writes)
        size_limit = (self.current_path.stat().st_size >= self.max_size)
        return (write_limit or size_limit)

    def get_file(self):
        """Return current HDF5 file, creating if necessary."""
        # Create new file if necessary
        if os.path.exists(str(self.current_path)):
            if self.check_file_limits():
                self.set_num += 1
                self.create_current_file()
        else:
            self.create_current_file()
        # Open current file
        comm = self.domain.distributor.comm_cart
        h5file = h5py.File(str(self.current_path), 'a', driver='mpio', comm=comm)
        self.file_write_num = h5file['/scales/write_number'].shape[0]
        return h5file

    @property
    def current_path(self):
        domain = self.domain
        comm = domain.distributor.comm_cart
        set_num = self.set_num
        file_name = '%s_s%i.hdf5' %(self.base_path.stem, set_num)
        return self.base_path.joinpath(file_name)

    def create_current_file(self):
        """Generate new HDF5 file in current_path."""
        self.file_write_num = 0
        comm = self.domain.distributor.comm_cart
        file = h5py.File(str(self.current_path), 'w-', driver='mpio', comm=comm)
        self.setup_file(file)

    def setup_file(self, file):

        domain = self.domain

        # Metadeta
        file.attrs['set_number'] = self.set_num
        file.attrs['handler_name'] = self.base_path.stem
        file.attrs['writes'] = self.file_write_num

        # Scales
        scale_group = file.create_group('scales')
        # Start time scales with shape=(0,) to chunk across writes
        scale_group.create_dataset(name='sim_time', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='timestep', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='world_time', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='wall_time', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='iteration', shape=(0,), maxshape=(None,), dtype=np.int)
        scale_group.create_dataset(name='write_number', shape=(0,), maxshape=(None,), dtype=np.int)
        const = scale_group.create_dataset(name='constant', data=np.array([0.], dtype=np.float64))
        for axis, basis in enumerate(domain.bases):
            coeff_name = basis.element_label + basis.name
            scale_group.create_dataset(name=coeff_name, data=basis.elements)
            scale_group.create_group(basis.name)

        # Tasks
        task_group =  file.create_group('tasks')
        for task_num, task in enumerate(self.tasks):
            layout = task['layout']
            scales = task['scales']
            gnc_shape, gnc_start, write_shape, write_start, write_count = self.get_write_stats(layout, scales, index=0)
            if np.prod(write_shape) <= 1:
                # Start with shape[0] = 0 to chunk across writes for scalars
                file_shape = (0,) + tuple(write_shape)
            else:
                # Start with shape[0] = 1 to chunk within writes
                file_shape = (1,) + tuple(write_shape)
            file_max = (None,) + tuple(write_shape)
            dset = task_group.create_dataset(name=task['name'], shape=file_shape, maxshape=file_max, dtype=layout.dtype)

            # Metadata and scales
            dset.attrs['task_number'] = task_num
            dset.attrs['grid_space'] = layout.grid_space
            dset.attrs['scales'] = scales

            # Time scales
            dset.dims[0].label = 't'
            for sn in ['sim_time', 'world_time', 'wall_time', 'timestep', 'iteration', 'write_number']:
                scale = scale_group[sn]
                dset.dims.create_scale(scale, sn)
                dset.dims[0].attach_scale(scale)

            # Spatial scales
            for axis, basis in enumerate(domain.bases):
                if layout.grid_space[axis]:
                    sn = basis.name
                    axscale = scales[axis]
                    if str(axscale) not in scale_group[sn]:
                        scale_group[sn].create_dataset(name=str(axscale), data=basis.grid(axscale))
                    lookup = '/'.join((sn, str(axscale)))
                else:
                    sn = lookup = basis.element_label + basis.name
                scale = scale_group[lookup]
                dset.dims.create_scale(scale, lookup)
                dset.dims[axis+1].label = sn
                dset.dims[axis+1].attach_scale(scale)

    def process(self, world_time, wall_time, sim_time, timestep, iteration, **kw):
        """Save task outputs to HDF5 file."""

        file = self.get_file()
        self.total_write_num += 1
        self.file_write_num  += 1
        file.attrs['writes'] = self.file_write_num
        index = self.file_write_num - 1

        # Update time scales
        sim_time_dset   = file['scales/sim_time']
        world_time_dset = file['scales/world_time']
        wall_time_dset  = file['scales/wall_time']
        timestep_dset   = file['scales/timestep']
        iteration_dset  = file['scales/iteration']
        write_num_dset  = file['scales/write_number']

        sim_time_dset.resize(index+1, axis=0)
        sim_time_dset[index] = sim_time
        world_time_dset.resize(index+1, axis=0)
        world_time_dset[index] = world_time
        wall_time_dset.resize(index+1, axis=0)
        wall_time_dset[index] = wall_time
        timestep_dset.resize(index+1, axis=0)
        timestep_dset[index] = timestep
        iteration_dset.resize(index+1, axis=0)
        iteration_dset[index] = iteration
        write_num_dset.resize(index+1, axis=0)
        write_num_dset[index] = self.total_write_num

        # Create task datasets
        for task_num, task in enumerate(self.tasks):
            out = task['out']
            out.set_scales(task['scales'], keep_data=True)
            out.require_layout(task['layout'])

            dset = file['tasks'][task['name']]
            dset.resize(index+1, axis=0)

            memory_space, file_space = self.get_hdf5_spaces(out.layout, task['scales'], index)
            dset.id.write(memory_space, file_space, out.data, dxpl=self._property_list)
            out.data=None
        file.close()

    def get_write_stats(self, layout, scales, index):
        """Determine write parameters for nonconstant subspace of a field."""

        # References
        gshape = layout.global_shape(scales)
        lshape = layout.local_shape(scales)
        start = layout.start(scales)
        first = (start == 0)

        # Build counts, taking just the first entry along constant axes
        write_count = lshape.copy()

        # Collectively writing global data
        global_nc_shape = gshape.copy()
        global_nc_start = start.copy()

        # Collectively writing global data
        write_shape = global_nc_shape
        write_start = global_nc_start

        return global_nc_shape, global_nc_start, write_shape, write_start, write_count

    def get_hdf5_spaces(self, layout, scales, index):
        """Create HDF5 space objects for writing nonconstant subspace of a field."""

        # References
        lshape = layout.local_shape(scales)
        start = layout.start(scales)
        gnc_shape, gnc_start, write_shape, write_start, write_count = self.get_write_stats(layout, scales, index)

        # Build HDF5 spaces
        memory_shape = tuple(lshape)
        memory_start = tuple(0 * start)
        memory_count = tuple(write_count)
        memory_space = h5py.h5s.create_simple(memory_shape)
        memory_space.select_hyperslab(memory_start, memory_count)

        file_shape = (index+1,) + tuple(write_shape)
        file_start = (index,) + tuple(write_start)
        file_count = (1,) + tuple(write_count)
        file_space = h5py.h5s.create_simple(file_shape)
        file_space.select_hyperslab(file_start, file_count)

        return memory_space, file_space

  
# Each node accumulates its own values
def add_t(self,solver,problem,dt):
    var=problem.variables

    u=None
    v=None
    w=None
    b=None
    s=None
    p=None

    a=dt*self.dy
    
    if 'u' in solver.state['u']: u=solver.state['u']['g']
    if 'v' in solver.state['v']: v=solver.state['v']['g']
    if 'w' in solver.state['w']: w=solver.state['w']['g']
    if 'b' in solver.state['b']: b=solver.state['b']['g']
    if 's' in solver.state['s']: s=solver.state['s']['g']
    if 'p' in solver.state['b']: p=solver.state['p']['g']

    if 'p'  in self: self.p  += a*np.sum(p,axis=1).reshape(self.sz)
    if 'b'  in self: self.b  += a*np.sum(b,axis=1).reshape(self.sz)
    if 's'  in self: self.s  += a*np.sum(s,axis=1).reshape(self.sz)
    if 'u'  in self: self.u  += a*np.sum(u,axis=1).reshape(self.sz)
    if 'v'  in self: self.v  += a*np.sum(v,axis=1).reshape(self.sz)
    if 'w'  in self: self.w  += a*np.sum(w,axis=1).reshape(self.sz)
    if 'bu' in self: self.bu += a*np.sum(b*u,axis=1).reshape(self.sz)
    if 'bv' in self: self.bv += a*np.sum(b*v,axis=1).reshape(self.sz)
    if 'bw' in self: self.bw += a*np.sum(b*w,axis=1).reshape(self.sz)
    if 'su' in self: self.su += a*np.sum(s*u,axis=1).reshape(self.sz)
    if 'sv' in self: self.sv += a*np.sum(s*v,axis=1).reshape(self.sz)
    if 'sw' in self: self.sw += a*np.sum(s*w,axis=1).reshape(self.sz)
    if 'uu' in self: self.uu += a*np.sum(u**2,axis=1).reshape(self.sz)
    if 'vv' in self: self.vv += a*np.sum(v**2,axis=1).reshape(self.sz)
    if 'ww' in self: self.ww += a*np.sum(w**2,axis=1).reshape(self.sz)
    if 'uv' in self: self.uv += a*np.sum(u*u,axis=1).reshape(self.sz)
    if 'uw' in self: self.uw += a*np.sum(u*v,axis=1).reshape(self.sz)
    if 'vw' in self: self.vw += a*np.sum(v*w,axis=1).reshape(self.sz)
    

def new(self,fn,t):
    self.t0=t
    
def save(self,fn,t):
    self.t1=t

    layout = domain.dist.grid_layout
    gshape = layout.global_shape(scales)
    lshape = layout.local_shape(scales)
    start  = layout.start(scales)

    # This 
    for i in range(gshape[0]):
        for j in range(gshape[2]):
            y[i,j]=comm.allreduce(y[i,j], op=MPI.SUM)
            

from setup_file in evaluator.py

def create_and_setup_file(self):
        """Generate new HDF5 file in current_path."""
        self.file_write_num = 0
        comm = self.domain.distributor.comm_cart
        file = h5py.File(str(self.current_path), 'w-', driver='mpio', comm=comm)
       domain = self.domain

        # Metadeta
        file.attrs['set_number'] = self.set_num
        file.attrs['handler_name'] = self.base_path.stem
        file.attrs['writes'] = self.file_write_num
        if not self.parallel:
            file.attrs['mpi_rank'] = domain.distributor.comm_cart.rank
            file.attrs['mpi_size'] = domain.distributor.comm_cart.size

        # Scales
        scale_group = file.create_group('scales')
        # Start time scales with shape=(0,) to chunk across writes
        scale_group.create_dataset(name='sim_time', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='timestep', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='world_time', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='wall_time', shape=(0,), maxshape=(None,), dtype=np.float64)
        scale_group.create_dataset(name='iteration', shape=(0,), maxshape=(None,), dtype=np.int)
        scale_group.create_dataset(name='write_number', shape=(0,), maxshape=(None,), dtype=np.int)
        const = scale_group.create_dataset(name='constant', data=np.array([0.], dtype=np.float64))
        for axis, basis in enumerate(domain.bases):
            coeff_name = basis.element_label + basis.name
            scale_group.create_dataset(name=coeff_name, data=basis.elements)
            scale_group.create_group(basis.name)

        # Tasks
        task_group =  file.create_group('tasks')
        for task_num, task in enumerate(self.tasks):
            layout = task['layout']
            constant = task['operator'].meta[:]['constant']
            scales = task['scales']
            gnc_shape, gnc_start, write_shape, write_start, write_count = self.get_write_stats(layout, scales, constant, index=0)
            if np.prod(write_shape) <= 1:
                # Start with shape[0] = 0 to chunk across writes for scalars
                file_shape = (0,) + tuple(write_shape)
            else:
                # Start with shape[0] = 1 to chunk within writes
                file_shape = (1,) + tuple(write_shape)
            file_max = (None,) + tuple(write_shape)
            dset = task_group.create_dataset(name=task['name'], shape=file_shape, maxshape=file_max, dtype=layout.dtype)
            if not self.parallel:
                dset.attrs['global_shape'] = gnc_shape
                dset.attrs['start'] = gnc_start
                dset.attrs['count'] = write_count

            # Metadata and scales
            dset.attrs['task_number'] = task_num
            dset.attrs['constant'] = constant
            dset.attrs['grid_space'] = layout.grid_space
            dset.attrs['scales'] = scales

            # Time scales
            dset.dims[0].label = 't'
            for sn in ['sim_time', 'world_time', 'wall_time', 'timestep', 'iteration', 'write_number']:
                scale = scale_group[sn]
                dset.dims.create_scale(scale, sn)
                dset.dims[0].attach_scale(scale)

            # Spatial scales
            for axis, basis in enumerate(domain.bases):
                if constant[axis]:
                    sn = lookup = 'constant'
                else:
                    if layout.grid_space[axis]:
                        sn = basis.name
                        axscale = scales[axis]
                        if str(axscale) not in scale_group[sn]:
                            scale_group[sn].create_dataset(name=str(axscale), data=basis.grid(axscale))
                        lookup = '/'.join((sn, str(axscale)))
                    else:
                        sn = lookup = basis.element_label + basis.name
                scale = scale_group[lookup]
                dset.dims.create_scale(scale, lookup)
                dset.dims[axis+1].label = sn
                dset.dims[axis+1].attach_scale(scale)


        
