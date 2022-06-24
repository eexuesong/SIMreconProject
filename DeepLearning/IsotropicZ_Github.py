from __future__ import print_function, unicode_literals, absolute_import, division
import time
import os
# this script was written by Yicong Wu for improving axial resolution in fluorescence microscopy on 06/17/2022
import numpy as np
import matplotlib.pyplot as plt
## matplotlib inline


##from tifffile import imread
import tifffile as tiff
from csbdeep.utils import axes_dict, plot_some, plot_history
from csbdeep.utils.tf import limit_gpu_memory
from csbdeep.io import load_training_data, save_tiff_imagej_compatible
from csbdeep.models import Config, CARE
from csbdeep.data import RawData, create_patches
from scipy.ndimage import zoom
from scipy.ndimage import rotate
from scipy.ndimage import gaussian_filter1d

#mode ='training' ## or set as 'validating' with ground_truth, 'testing' without ground_truth
mode='prediction' ## or set as 'validating' with ground_truth, 'testing' without ground_truth

## folders for training
training_dir = 'Z:\\Xuesong\\3D_SIM_DL\\2022_02_02_001_Jurkat_EMTB_EGFP_high_low_SNR(Nikon)\\Wiener_reconstruction(highSNR)\\50nm'
input_folders = ['XZ_1DSIM_Input']
truth_folder  = 'XZ_1DSIM_GT'

## folders for the model to call for recovery
model_dir ='Z:\\Xuesong\\3D_SIM_DL\\all models\\Lamp1_Model'
model_name = 'CARE_XZ_6degree_Model_Lamp1_50nm'

## folders for prediction
prediction_dir = 'Z:\\Xuesong\\3D_SIM\\2022_03_15_001_Lamp1_GFP_Lysotracker_Red_RT(Nikon)\\lamp_Wiener_reconstruction_001_denoisingStep2'

if mode == 'prediction':
    output_folder = prediction_dir + '_CARE_6Degree_DL'
    try:
        DL_path = output_folder
        if not os.path.exists(DL_path):
            os.makedirs(DL_path)
    except OSError:
        print ("Creation of the directory %s failed" % DL_path)
    else:
        print ("Successfully created the directory %s " % DL_path)

# change the image volume to 50 nm pixel size in x, y, z
zoom_x = 40.88/50
zoom_y = zoom_x
zoom_z = 125/50

# set preproces_flag
preproces_flag = 1

# set the blurring kernel
sigma_x = 2.6
sigma_y = 1.3  #1.3;  0.6

# set zero pad or not for rotation: 1 yes, 0 no
pad_flag = 1

# set saving intermediate files or not: 1 yes, 0 no
save_flag = 0

# set upsampling/downsampling zoom flag after blurring
resampling_flag = 1

# set training size

if mode == 'training':
    gt_labels = [f for f in os.listdir(training_dir + '\\' + truth_folder) if os.path.isfile(os.path.join(training_dir + '\\' + truth_folder, f))]
    gt_File = training_dir +'\\'+ truth_folder + '\\'+ gt_labels[0]
    gt_data = tiff.imread(gt_File)
    (nz_gt, ny_gt, nx_gt) = gt_data.shape
    print('GT training data size:', (nz_gt, ny_gt, nx_gt))

    # patch size must be 4 x
    if nz_gt <= 64:
        patch_sizeZ = nz_gt
    else:
        patch_sizeZ = 64

    if ny_gt <= 64:
        patch_sizeY = ny_gt
    else:
        patch_sizeY = 64

    if nx_gt <= 64:
        patch_sizeX = nx_gt
    else:
        patch_sizeX = 64

    patch_size = (patch_sizeZ, patch_sizeY, patch_sizeX) # z, y, x (40,64,64)
    patch_number = int(nz_gt*ny_gt*nx_gt/(patch_sizeZ*patch_sizeY*patch_sizeX)) # 300 for z need to calculate based on the volume size
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

def training():
    raw_data = RawData.from_folder(
    basepath    = training_dir,
    source_dirs = input_folders,
    target_dir  = truth_folder,
    axes        = 'ZYX',
    )

    X, Y, XY_axes = create_patches(
    raw_data            = raw_data,
    patch_size          = patch_size,
    n_patches_per_image = patch_number,
    save_file           = training_dir + '\\' + model_name + '.npz'
    )

    assert X.shape == Y.shape
    print("shape of X,Y =", X.shape)
    print("axes  of X,Y =", XY_axes)

    (X,Y), (X_val,Y_val), axes = load_training_data(training_dir + '\\' + model_name + '.npz', validation_split=0.1, verbose=True)
    c = axes_dict(axes)['C']
    n_channel_in, n_channel_out = X.shape[c], Y.shape[c]

    config = Config(axes, n_channel_in, n_channel_out, train_epochs = 100, train_steps_per_epoch=100)
    print(config)
    vars(config)
    model = CARE(config, model_name, basedir=training_dir)
    history = model.train(X,Y, validation_data=(X_val,Y_val))
    print(sorted(list(history.history.keys())))
    plt.figure(figsize=(16,5))
    plot_history(history,['loss','val_loss'],['mse','val_mse','mae','val_mae']);

def prediction():
    axes = 'ZYX'
    #axes = 'YX'
    decon_model = CARE(config=None, name=model_name, basedir=model_dir)

    input_labels = [f for f in os.listdir(prediction_dir) if os.path.isfile(os.path.join(prediction_dir, f))]

    print (input_labels)
    maxlen = len(input_labels)

    angle = np.linspace(-90, 90, 7)

    for i in range(0,maxlen):
        print('processing.... ' + input_labels[i])
        Prediction_File = prediction_dir + '\\'+ input_labels[i]
        input_data = tiff.imread(Prediction_File)

        if preproces_flag == 1:
            input_data = zoom(input_data, (zoom_z, zoom_y, zoom_x), mode='reflect')
            input_data = np.swapaxes(input_data,0,2)
            input_data = np.rot90(input_data, k=1, axes=(1, 2))
            input_data = np.flip(input_data, axis=1)
            tiff.imsave(output_folder + '\\50nm_' + input_labels[i], input_data)
            input_data = gaussian_filter1d(input_data, sigma_x, axis=2)
            input_data = gaussian_filter1d(input_data, sigma_y, axis=1)
            if resampling_flag == 1:
                input_data = zoom(input_data, (1, 1, 1/zoom_z), mode='reflect')
                input_data = zoom(input_data, (1, 1, zoom_z), mode='reflect')

            if save_flag == 1:
                tiff.imsave(output_folder + '\\blur_50nm_' + input_labels[i], input_data)
        input_data[input_data<=0]=0
        (nz, ny, nx) = input_data.shape
        results = np.zeros((len(angle) - 1, nz, ny, nx))
        print(results.shape)

        tile_sizeX = ny
        overlap_size = 10

        tile_x = nx/tile_sizeX
        if tile_x == 1:
            tile_x = int(tile_x)
        else:
            tile_x = int(tile_x) + 1

        start_x = 0
        for xt in range(tile_x):
            print('number of tiles:', tile_x, '... processing tile #', xt+1)
            start_x = xt * tile_sizeX
            end_x = min(start_x + tile_sizeX, nx)
            tile_data = input_data[:, :, max(start_x-overlap_size,0): min(end_x+overlap_size,nx)]
            (nzt, nyt, nxt) = tile_data.shape
            print('tile_data:', tile_data.shape)
            nx1 = max(int(np.sqrt(nxt * nxt + nyt * nyt)),int(nyt*1.414))
            pad_y = round((nx1 - nyt) / 2)
            pad_x = round((nx1 - nxt) / 2)
            tile_data_pad = np.pad(tile_data, ((0, 0), (pad_y, nx1-nyt-pad_y), (pad_x, nx1-nxt-pad_x)), mode='reflect')
            [nz2, ny2, nx2] = tile_data_pad.shape
            print('tile_data_pad:', tile_data_pad.shape)
            # tiff.imsave(output_folder + '\\pad_' + input_labels[i], input_data_pad)
            for k in range(len(angle)-1):
           # rotation
                if angle[k] == 0:
                    tile_results = decon_model.predict(tile_data, axes)
                    print(str(int(angle[k])), tile_results.shape)
                elif angle[k] == -90:
                    output_data = np.rot90(tile_data, k=1, axes=(1, 2))
                    output_data = decon_model.predict(output_data, axes)
                    tile_results = np.rot90(output_data, k=1, axes=(2, 1))
                    print(str(int(angle[k])), tile_results.shape)
                else:
                    output_data = rotate(tile_data_pad, -1 * angle[k], axes=(1, 2), reshape=False)
                    output_data = decon_model.predict(output_data, axes, n_tiles=(1, 1, 1))
                    output_data = rotate(output_data, angle[k], axes=(1,2), reshape=False)
                    tile_results = output_data[:, pad_y:pad_y + nyt, pad_x:pad_x + nxt]
                    print(str(int(angle[k])), tile_results.shape)
                nx3 = end_x - start_x
                if xt == 0:
                    new_start_x = 0
                else:
                    new_start_x = overlap_size
                results[k,:,:,start_x:end_x] = tile_results[:,:,new_start_x:nx3+new_start_x]

           # normalize and save
           #  inmin, inmax = final_result.min(), final_result.max()
           #  final_result = (final_result - inmin) / (inmax - inmin) * 65535
        if save_flag == 1:
            for k in range(len(angle) - 1):
                tiff.imsave(output_folder + '\\DL_' + str(int(angle[k])) + '_' + input_labels[i], results[k, :, :, :].astype('single'))
        results[results <= 0] = 0.000001

        print('DL results sahpe:', results.shape)
        mean_dir = np.mean(results,axis=(1,2,3))
        mean_max = np.amax(mean_dir)
        results_fft = np.zeros_like(results,dtype=complex)
        for k in range(len(angle)-1):
            results[k, :, :, :] = results[k, :, :, :] * mean_max / mean_dir[k]
            results_fft[k,:,:,:] = np.fft.fftn(results[k,:,:,:])
        results_fft_abs = np.absolute(results_fft)

        max_index = np.argmax(results_fft_abs, axis=0)
        Wiener_fft_max = np.zeros_like(results_fft[0,:,:,:])

        for k in range(len(angle)-1):
            max_fft_abs = np.zeros_like(input_data)
            max_fft_abs[max_index==k] = 1
            Wiener_fft_max = Wiener_fft_max + results_fft[k,:,:,:] * max_fft_abs
        Wiener_max = np.fft.ifftn(Wiener_fft_max)
        Wiener_max_abs = np.absolute(Wiener_max).astype('single')

        print(Wiener_max_abs.shape)
        tiff.imsave(output_folder + '\\max_fft_abs_' + input_labels[i], Wiener_max_abs)

        if save_flag == 1:
            Wiener_max_real = np.real(Wiener_max).astype('single')
            tiff.imsave(output_folder + '\\max_fft_real_' + input_labels[i], Wiener_max_real)


def main():
    print('CARE') 
    if mode=='training':
        training()
    elif mode=='validation':
        validation()
    elif mode=='prediction':
        prediction()
    else: print('wrong mode!')

main()

