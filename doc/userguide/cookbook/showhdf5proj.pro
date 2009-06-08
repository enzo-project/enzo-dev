 pro showhdf5proj, file_name, outfile_name, number

file_id = H5F_OPEN(file_name)


if (number eq 4) then dataset_id = H5D_OPEN(file_id,'/projected_gas_density')
if (number eq 7) then dataset_id = H5D_OPEN(file_id,'/projected_x-ray_luminosity_div1e23')
if (number eq 3) then dataset_id = H5D_OPEN(file_id,'/projected_dm_density')
if (number eq 8) then dataset_id = H5D_OPEN(file_id,'/projected_x-ray_weighted_temperature')
if (number eq 5) then dataset_id = H5D_OPEN(file_id,'/projected_level')
if (number eq 2) then dataset_id = H5D_OPEN(file_id,'/SZ_y_effect')
if (number eq 0) then dataset_id = H5D_OPEN(file_id,'/DT_over_T_Doppler_effect')
if (number eq 1) then dataset_id = H5D_OPEN(file_id,'/Metallicity')
if (number eq 6) then dataset_id = H5D_OPEN(file_id,'/projected_star_density')

print, 'test ', file_name, number
print, 'dataset id ', dataset_id

data = H5D_READ(dataset_id)

nummem = H5G_GET_NMEMBERS(file_id,'/')

print, 'there are ',nummem,' members'
membername = H5G_GET_MEMBER_NAME(file_id,'/',number)
print, 'we are opening ',membername

dataspace_id = H5D_GET_SPACE(dataset_id)


dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

print, 'dimensions ', dimensions

H5S_CLOSE, dataspace_id
H5D_CLOSE, dataset_id
H5F_CLOSE, file_id


;WINDOW, XSIZE=dimensions[0], YSIZE=dimensions[1]

set_plot,'z'
device,set_resolution=[dimensions[0],dimensions[1]]

data_dim = size(data)
print, 'image info: ',data_dim
;window

if (number eq 0) then data2 = data
if (number eq 1) then data2 = ALOG10(data)
if (number eq 2) then data2 = data
if (number eq 3) then data2 = ALOG10(data+1.0e12)
if (number eq 4) then data2 = ALOG10(data)
;if (number eq 3) then data2 = data
if (number eq 5) then data2 = data
if (number eq 6) then data2 = ALOG10(data+1.0e9) 
if (number eq 7) then data2 = ALOG10(data)
if (number eq 8) then data2 = ALOG10(data)


;min1 = MIN(data2)
;max1 = MAX(data2)

;print,  'min, max =', min1,max1

;surface,data2


TVSCL, data2

image=TVRD()
image3D = bytarr(3,dimensions[0],dimensions[1]) 
tvlct,r,g,b,/get

image3d[0,*,*] = r[image]
image3d[1,*,*] = g[image]
image3d[2,*,*] = b[image]

write_jpeg,outfile_name,image3d,true=1,quality=100


END 
