def gather_field(field):

  field.require_layout('c')
  slices = field.layout.slices((1,1,1))

  array = np.zeros(field.layout.global_shape((1,1,1)))

  data = domain.distributor.comm_world.gather(field['c'], root=0)
  sliceses = domain.distributor.comm_world.gather(slices, root=0)

  if domain.distributor.rank==0:
    for (slices,datum) in zip(sliceses,data):
      array[np.array(slices)] = datum
  else:
    array = None

  array = domain.distributor.comm_world.scatter([array]*domain.distributor.size, root=0)

  return array
