import numpy as np
from sys import argv
from os import listdir
from os.path import isfile, join
import os, json
import argparse

from utils import read_config, get_data, get_patch

def main(argv):

    parser = argparse.ArgumentParser(description='Makes training data set for EM vs track separation')
    parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
    args = parser.parse_args()

    config = read_config(args.config)

    print '#'*50,'\nPrepare data for CNN'
    INPUT_DIR = config['prepare_data_em_track']['input_dir']
    OUTPUT_DIR = config['prepare_data_em_track']['output_dir']
    PATCH_SIZE_W = config['prepare_data_em_track']['patch_size_w']
    PATCH_SIZE_D = config['prepare_data_em_track']['patch_size_d']
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    doing_nue = config['prepare_data_em_track']['doing_nue']                       # set to true for nu_e events (will skip more showers)
    selected_view_idx = config['prepare_data_em_track']['selected_view_idx']       # set the view id
    patch_fraction = config['prepare_data_em_track']['patch_fraction']             # percent of used patches
    empty_fraction = config['prepare_data_em_track']['empty_fraction']             # percent of "empty background" patches
    clean_track_fraction = config['prepare_data_em_track']['clean_track_fraction'] # percent of selected patches, where only a clean track is present
    crop_event = config['prepare_data_em_track']['crop_event']                     # use true only if no crop on LArSoft level and not a noise dump

    print 'Using', patch_fraction, '% of data from view', selected_view_idx
    if doing_nue: print 'Neutrino mode, will skip more showers.'

    max_capacity = 1700000
    db = np.zeros((max_capacity, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
    db_y = np.zeros((max_capacity, 3), dtype=np.int32)

    patch_area = PATCH_SIZE_W * PATCH_SIZE_D

    cnt_ind = 0
    cnt_trk = 0
    cnt_sh = 0
    cnt_void = 0

    fcount = 0

    files = [f for f in os.listdir(INPUT_DIR) if '.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        finfo = fname.split('_')
        evt_no = finfo[2]
        tpc_idx = int(finfo[8])
        view_idx = int(finfo[10])

        if view_idx != selected_view_idx: continue
        fcount += 1

        print 'Process file', fcount, fname, 'EVT', evt_no

        # get clipped data, margin depends on patch size in drift direction
        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname, PATCH_SIZE_D/2 + 2, crop_event)
        if raw is None:
            print 'Skip empty event...'
            continue

        print 'Tracks', np.sum(tracks), 'showers', np.sum(showers)

        sel_trk = 0
        sel_sh = 0
        for i in range(raw.shape[0]):
            for j in range(raw.shape[1]):
                # randomly skip fraction of patches
                if np.random.randint(10000) > int(100*patch_fraction): continue
                
                x_start = np.max([0, i - PATCH_SIZE_W/2])
                x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE_W])

                y_start = np.max([0, j - PATCH_SIZE_D/2])
                y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE_D])

                if x_stop - x_start != PATCH_SIZE_W or y_stop - y_start != PATCH_SIZE_D:
                    continue

                track_pixels = np.count_nonzero(tracks[x_start:x_stop, y_start:y_stop])
                shower_pixels = np.count_nonzero(showers[x_start:x_stop, y_start:y_stop])

                target = np.zeros(3)
                if tracks[i,j] == 1:
                    # skip some fraction (1/2) of almost-track-only patches
                    if shower_pixels < 4 and np.random.randint(10000) < int(100*clean_track_fraction): continue
                    else:
                        target[0] = 1
                        cnt_trk += 1
                        sel_trk += 1
                elif showers[i,j] == 1:
                    # for nu_e events (lots of showers) skip some fraction of shower patches
                    if doing_nue:
                        if np.random.randint(100) < 40: continue  # skip 40% of any shower
                        if shower_pixels > 0.05*patch_area and np.random.randint(100) < 50: continue
                        if shower_pixels > 0.20*patch_area and np.random.randint(100) < 90: continue
                    target[1] = 1
                    cnt_sh += 1
                    sel_sh += 1
                else:
                    # use small fraction (1/100) of empty-but-close-to-something patches
                    if np.random.randint(10000) < int(100*empty_fraction):
                        nclose = np.count_nonzero(showers[i-2:i+3, j-2:j+3])
                        nclose += np.count_nonzero(tracks[i-2:i+3, j-2:j+3])
                        if nclose == 0:
                            npix = shower_pixels + track_pixels
                            if npix > 6:
                                target[2] = 1
                                cnt_void += 1
                            else: continue # completely empty patch
                        else: continue # too close to trk/shower
                    else: continue # not selected randomly

                if cnt_ind < max_capacity:
                    db[cnt_ind] = get_patch(raw, i, j, PATCH_SIZE_W, PATCH_SIZE_D)
                    db_y[cnt_ind] = target
                    cnt_ind += 1
                else: break

        print 'Tracks', sel_trk, 'showers', sel_sh, 'selected'

    print 'Added', cnt_ind, 'Tracks:', cnt_trk, 'showers:', cnt_sh, 'empty:', cnt_void

    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_x', db[:cnt_ind])
    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_y', db_y[:cnt_ind])

if __name__ == "__main__":
    main(argv)

