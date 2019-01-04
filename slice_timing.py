from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

import nibabel as nib
img = nib.load('an_example_4d.nii')
data = img.get_data()

data.shape

vol0 = data[..., 0]

plt.rcParams['image.cmap'] = 'gray'  
plt.rcParams['image.interpolation'] = 'nearest'

plt.imshow(vol0[31, :, :].T, origin='bottom left')  
plt.title('Sagittal section through first volume') 
plt.xlabel('x axis')  
plt.ylabel('z axis')  

TR = 2.0
n_z_slices = 16
time_for_single_slice = TR / n_z_slices
time_for_single_slice

time_for_slice_0 = 0
time_for_slice_1 = time_for_single_slice * 8
time_for_slice_1

plt.imshow(vol0[:, :, 0])  
plt.title('Vol 0, z slice 0')  

vox_x = 14  
vox_y = 22  

time_course_slice_0 = data[vox_x, vox_y, 0, :]
time_course_slice_1 = data[vox_x, vox_y, 1, :]

vol_nos = np.arange(data.shape[-1])
vol_onset_times = vol_nos * TR
times_slice_0 = vol_onset_times

times_slice_1 = vol_onset_times + TR / 2.

plt.plot(times_slice_0, time_course_slice_0, 'b:+',
    label='slice 0 time course')
plt.plot(times_slice_1, time_course_slice_1, 'r:+',
    label='slice 1 time course')
plt.legend()
plt.title('Time courses for slice 0, slice 1')
plt.xlabel('time (seconds)')  

plt.plot(times_slice_0[:10], time_course_slice_0[:10], 'b:+',
    label='slice 0 time course')
plt.plot(times_slice_1[:10], time_course_slice_1[:10], 'r:+',
    label='slice 1 time course')
plt.legend()
plt.title('First 10 values for slice 0, slice 1')
plt.xlabel('time (seconds)')  

plt.plot(times_slice_0[:10], time_course_slice_0[:10], 'b:+')
plt.plot(times_slice_1[:10], time_course_slice_1[:10], 'r:+')
plt.title('First 10 values for slice 0, slice 1')
plt.xlabel('time (seconds)')  
min_y, max_y = plt.ylim()
for i in range(1, 10):
    t = times_slice_0[i]
    plt.plot([t, t], [min_y, max_y], 'k:')

plt.plot(times_slice_0[:10], time_course_slice_0[:10], 'b:+')
plt.plot(times_slice_1[:10], time_course_slice_1[:10], 'r:+')
plt.title('First 10 values for slice 0, slice 1')
plt.xlabel('time (seconds)')
min_y, max_y = plt.ylim()
for i in range(1, 10):
    t = times_slice_0[i]
    plt.plot([t, t], [min_y, max_y], 'k:')
    x = t
    x0 = times_slice_1[i-1]
    x1 = times_slice_1[i]
    y0 = time_course_slice_1[i-1]
    y1 = time_course_slice_1[i]

    y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    plt.plot(x, y, 'kx')

from scipy.interpolate import InterpolatedUnivariateSpline as Interp

lin_interper = Interp(times_slice_1, time_course_slice_1, k=1)
type(lin_interper)

interped_vals = lin_interper(times_slice_0)

plt.plot(times_slice_0[:10], time_course_slice_0[:10], 'b:+')
plt.plot(times_slice_1[:10], time_course_slice_1[:10], 'r:+')
plt.plot(times_slice_0[:10], interped_vals[:10], 'kx')
plt.title('Using the scipy interpolation object')

plt.plot(times_slice_0, interped_vals, 'r:+',
    label='interpolated slice 1 time course')
plt.plot(times_slice_0, time_course_slice_0, 'b:+',
    label='slice 0 time course')
plt.legend()
plt.title('Slice 1 time course interpolated to slice 0 times')
plt.xlabel('time (seconds)')  
