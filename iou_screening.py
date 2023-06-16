import os, sys
import pandas as pd
from tqdm import tqdm
from iou import iou_compute
import numpy as np


def main_process(path_true, path_pred, outpath, i, parallel=1, plt_scatter=False):
    j = int(i) + 10
    print(f'{i}-{j} begin!')
    pdb_ids = os.listdir(path_pred)
    atom_interact_iou_list, aa_interact_iou_list, key_aa_iou_list, screening_iou_list = [], [], [], []
    for pdb_id in tqdm(pdb_ids[int(i):j]):
        # ban_pdbids = ['1bzc', '1gpn', '5tmn']
        # if pdb_id in ban_pdbids: continue
        true_prot_lig_file = [os.path.join(path_true, pdb_id, f'{pdb_id}_protein.pdb'),
                              os.path.join(path_true, pdb_id, f'{pdb_id}_ligand.sdf')]
        if not os.path.exists(true_prot_lig_file[0]) or not os.path.exists(true_prot_lig_file[1]):
            print(f'{pdb_id}晶体文件不存在。')
            continue
        pair_data_name, pred_prot_lig_files = [], []
        for dir in os.listdir(os.path.join(path_pred, pdb_id)):
            # ligand_data = []
            for pose in os.listdir(os.path.join(path_pred, pdb_id, dir)):
                # if 'complex' in pose or pose.endswith('pdb'): continue
                if not pose.endswith('.sdf'): continue
                pair_data_name.append(pose[:-4])
                pred_prot_lig_files.append([None, os.path.join(path_pred, pdb_id, dir, pose)])
        pair_data = (true_prot_lig_file, pred_prot_lig_files)
        print(dir, '数量为:', len(pred_prot_lig_files))
        assert pair_data[0] != [] and pair_data[1][1] != []
        # pair_data = pair_data[:5]
        atom_interact_iou, aa_interact_iou, key_aa_iou, screening_iou = iou_compute(pair_data, parallel=parallel,
                                                                                         plt_scatter=plt_scatter)  # 目前merge_complex只支持pymol
        atom_interact_iou_list.extend(atom_interact_iou)
        aa_interact_iou_list.extend(aa_interact_iou)
        key_aa_iou_list.extend(key_aa_iou)
        screening_iou_list.extend(screening_iou)
        dic1 = {'#code_ligand_num': pair_data_name, 'score': aa_interact_iou}
        dic2 = {'#code_ligand_num': pair_data_name, 'score': key_aa_iou}
        dic3 = {'#code_ligand_num': pair_data_name, 'score': screening_iou}
        df1, df2, df3 = pd.DataFrame(dic1), pd.DataFrame(dic2), pd.DataFrame(dic3)
        os.system(
            f'mkdir -p {outpath}/1_aa_interaction_iou && mkdir -p {outpath}/2_key_aa_iou && mkdir -p {outpath}/3_screening_iou')
        df1.to_csv(f'{outpath}/1_aa_interaction_iou/{pdb_id}_score.dat', sep='\t', index=False)
        df2.to_csv(f'{outpath}/2_key_aa_iou/{pdb_id}_score.dat', sep='\t', index=False)
        df3.to_csv(f'{outpath}/3_screening_iou/{pdb_id}_score.dat', sep='\t', index=False)
    print('原子粒度IOU均值: ', np.mean(atom_interact_iou_list))
    print('残基粒度IOU均值: ', np.mean(aa_interact_iou_list))
    print('关键残基IOU均值: ', np.mean(key_aa_iou_list))
    print('关键残基复现比均值: ', np.mean(screening_iou_list))
    print(f'{i} - {i+9} have Done!')


if __name__ == '__main__':
    path_true = '/home/chentong/data/coreset'
    path_pred = '/home/chentong/data/decoys_screening'
    outpath = '/home/chentong/data/iou_screen'
    # path_true = 'example_data/true'
    # path_pred = 'example_data/pred'
    # outpath = 'example_data/ret'
    # path_true = '/mnt/d/interaction_iou/example_data/true'
    # path_pred = '/mnt/d/interaction_iou/example_data/pred'
    # outpath = '/mnt/d/data/CASF-2016/CASF-2016/000_iou'
    i = sys.argv[1] # 10个pdbid为一组，在一台服务器上跑    # i = 0, 10, 20, 30, 40, 50
    print(f'i的值为{i}')
    main_process(path_true, path_pred, outpath, i, parallel=10, plt_scatter=False)

