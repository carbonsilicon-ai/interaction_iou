import os
from tqdm import tqdm
import pandas as pd
from iou import iou_compute
import numpy as np


def main_process(path_true, path_pred, outpath, parallel=1, plt_scatter=False):
    pdb_ids = os.listdir(path_pred)
    atom_interact_iou_list, aa_interact_iou_list, key_aa_iou_list, screening_iou_list = [], [], [], []
    for pdb_id in tqdm(pdb_ids):
        # ban_pdbids = ['1bzc', '1gpn', '5tmn']
        # if pdb_id in ban_pdbids: continue
        true_prot_lig_file = [os.path.join(path_true, pdb_id, f'{pdb_id}_protein.pdb'),
                              os.path.join(path_true, pdb_id, f'{pdb_id}_ligand.sdf')]
        if not os.path.exists(true_prot_lig_file[0]) or not os.path.exists(true_prot_lig_file[1]):
            print(f'{pdb_id}晶体文件不存在。')
            continue
        pair_data_name, pred_prot_lig_files = [], []
        for dir in os.listdir(os.path.join(path_pred, pdb_id)):
            for pose in os.listdir(os.path.join(path_pred, pdb_id, dir)):
                # if 'complex' in pose or pose.endswith('pdb'): continue
                if not pose.endswith('.sdf'): continue
                pair_data_name.append(pose[:-4])
                pred_prot_lig_files.append([None, os.path.join(path_pred, pdb_id, dir, pose)])
        pair_data = (true_prot_lig_file, pred_prot_lig_files)
        print(dir, ':', len(pred_prot_lig_files))
        assert pair_data[0] != [] and pair_data[1][1] != []
        # pair_data = pair_data[:5]
        atom_interact_iou, aa_interact_iou, key_aa_iou, screening_iou = iou_compute(pair_data, parallel=parallel,
                                                                                         plt_scatter=plt_scatter)  # 目前merge_complex只支持pymol
        atom_interact_iou_list.extend(atom_interact_iou)
        aa_interact_iou_list.extend(aa_interact_iou)
        key_aa_iou_list.extend(key_aa_iou)
        screening_iou_list.extend(screening_iou)
    print('原子粒度IOU均值: ', np.mean(atom_interact_iou_list))
    print('残基粒度IOU均值: ', np.mean(aa_interact_iou_list))
    print('关键残基IOU均值: ', np.mean(key_aa_iou_list))
    print('关键残基复现比均值: ', np.mean(screening_iou_list))


if __name__ == '__main__':
    path_true = 'example_data/true'
    path_pred = 'example_data/pred'
    outpath = 'example_data/ret'
    main_process(path_true, path_pred, outpath, parallel=1, plt_scatter=False)
    print('Done!')

