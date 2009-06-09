filename = "GridDensity"
slice = 1024
resolution = 2048
output = "rho_slice.jpg"

fid = H5F_OPEN(filename)
data_id = H5D_OPEN(fid, "/GridDensity")
dspace_id = H5D_GET_SPACE(data_id)

count = [1,resolution,1,resolution]
H5S_SELECT_HYPERSLAB, dspace_id, [0,0,slice,0], count, /RESET
mspace_id = H5S_CREATE_SIMPLE(count)
data = H5D_READ(data_id, FILE_SPACE=dspace_id, MEMORY_SPACE=mspace_id)

H5S_CLOSE, mspace_id
H5S_CLOSE, dspace_id
H5D_CLOSE, data_id
H5F_CLOSE, fid

data = REFORM(data)

END
