import io
import os
# import prody
import time

from rdkit import Chem
import numpy as np
import pandas as pd
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport
import matplotlib.pyplot as plt
from pymol import cmd
from multiprocessing import Pool
from plip.basic import config


def retrieve_interactions(complex_pdb_str, pdb_file, merge_complex='pymol'):
    config.NOFIXFILE = True
    protlig = PDBComplex()
    pdb_dir = '_'.join(pdb_file.split('/'))
    path_tmp = os.path.join('/dev/shm', pdb_dir)
    os.system(f'mkdir -p {path_tmp}')
    protlig.output_path = path_tmp
    protlig.load_pdb(complex_pdb_str, as_string=True)
    if merge_complex == 'prody':
        for ligand in protlig.ligands:
            if ligand.longname != 'UNL': continue
            protlig.characterize_complex(ligand)
    else:
        for ligand in protlig.ligands:
            if ligand.longname != 'FFF': continue
            protlig.characterize_complex(ligand)

    sites = {}
    for key, site in sorted(protlig.interaction_sets.items()):
        binding_site = BindingSiteReport(site)  # collect data about interactions
        # tuples of *_features and *_info will be converted to pandas DataFrame
        keys = ("hydrophobic", "hbond", "waterbridge", "saltbridge", "pistacking", "pication", "halogen", "metal")
        interactions = {
            k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info")
            for k in keys
        }
        sites[key] = interactions
    os.system(f'rm -fr {protlig.output_path}')
    return sites


def create_df_from_binding_site(selected_site_interactions, interaction_type):
    df = pd.DataFrame.from_records(
        selected_site_interactions[interaction_type][1:],
        columns=selected_site_interactions[interaction_type][0],
    )
    if 'PROT_IDX_LIST' in df.columns:
        df['PROT_IDX_LIST'] = df['PROT_IDX_LIST'].apply(lambda x: ','.join(sorted(x.split(','))))
    if 'LIG_IDX_LIST' in df.columns:
        df['LIG_IDX_LIST'] = df['LIG_IDX_LIST'].apply(lambda x: ','.join(sorted(x.split(','))))
    return df


def interaction_extract(interactions_by_site):
    def apply_func(data):
        if type(data[0]) is str:
            if data[0].split(',')[0] > data[1].split(',')[0]:
                return data[1], data[0]
            else:
                return data[0], data[1]
        else:
            if data[0] > data[1]:
                return data[1], data[0]
            else:
                return data[0], data[1]

    interact_types = {"hydrophobic": ['PROTCARBONIDX', 'LIGCARBONIDX', 1],  # 111 通过名字区分蛋白和配体
                      "hbond": ['DONORIDX', 'ACCEPTORIDX', 'PROTISDON', 2],  # 111 通过PROTISDON来区分
                      "waterbridge": ['DONOR_IDX', 'ACCEPTOR_IDX', 'PROTISDON', 3],  # 111 通过PROTISDON来区分
                      "saltbridge": ['PROT_IDX_LIST', 'LIG_IDX_LIST', 4],  # 通过名字区分
                      "pistacking": ['PROT_IDX_LIST', 'LIG_IDX_LIST', 5],  # 通过名字区分
                      "pication": ['PROT_IDX_LIST', 'LIG_IDX_LIST', 6],  # 通过名字区分
                      "halogen": ['DON_IDX', 'ACC_IDX', 7],  # 111 DON_IDX即配体为卤键给体, ACC_IDX即蛋白为卤键受体
                      "metal": ['METAL_IDX', 'TARGET_IDX', 8]  # 111 METAL_IDX为蛋白，
                      }
    key_aa, aa_interaction, atom_interaction = [], [], []
    for selected_site, interact_ret in interactions_by_site.items():
        for interact_type, colnames in interact_types.items():
            df = create_df_from_binding_site(interact_ret, interaction_type=interact_type)
            if df.empty: continue
            if interact_type in ['hydrophobic', 'hbond', 'waterbridge', 'halogen', 'metal']:
                # hydro_pairs = zip(df['PROTCARBONIDX'].tolist(), df['LIGCARBONIDX'].tolist())
                # dics.extend(['_'.join((interact_type, str(lig), str(prot))) for prot, lig in hydro_pairs])
                new_data = df[colnames[:2]].apply(apply_func, axis=1).tolist()
            else:  # ['saltbridge', 'pistacking', 'pication']
                new_data = df[['PROT_IDX_LIST', 'LIG_IDX_LIST']].apply(apply_func, axis=1).tolist()
            atom_interaction.extend(['_'.join((interact_type, str(tmp[0]), str(tmp[1]))) for tmp in new_data])
            cur_aa = list(zip(df['RESNR'].tolist(), df['RESCHAIN'].tolist()))
            aa_interaction.extend(['_'.join((interact_type, str(aa[0]), str(aa[1]))) for aa in cur_aa])
            key_aa.extend(['_'.join((str(aa[0]), str(aa[1]))) for aa in cur_aa])
    return atom_interaction, aa_interaction, list(set(key_aa))


def interact_process(complex_pdb_str, pdb_file, merge_complex):
    interactions_by_site = retrieve_interactions(complex_pdb_str, pdb_file, merge_complex)
    atom_interaction, aa_interaction, key_aa = interaction_extract(interactions_by_site)
    return atom_interaction, aa_interaction, key_aa


def iou_metrics(true, pred):
    iou_ret = len(set(true) & (set(pred))) * 1.0 / (len(set(true) | set(pred)) + 1e-8)
    return iou_ret


def iou_virtual_screen_metrics(true, pred):
    iou_virtual_ret = len(set(true) & (set(pred))) * 1.0 / (len(set(true)) + 1e-8)
    return iou_virtual_ret


def build_complex_core_pymol(pair):
    protein_file, ligand_file = pair
    cmd.reinitialize()
    cmd.load(ligand_file, 'ligand')
    cmd.alter('ligand', 'resn="FFF"')
    cmd.remove('hydrogens')
    # cmd.sort('ligand')
    cmd.load(protein_file, "prot")
    cmd.set('pdb_use_ter_records', 0)
    cmd.create("complex", "prot, ligand")
    cmd.sort("ligand, prot")
    pdb_str = cmd.get_str('pdb', 'complex')
    # cmd.save(f, 'complex')
    return pdb_str


class MultiProcessing2(object):
    def __init__(self, merge_complex='pymol'):
        self.merge_complex = merge_complex

    def get_interaction(self, pair):
        complex_pdb_str = build_complex_core_pymol(pair)
        interactions = interact_process(complex_pdb_str, pair[1], self.merge_complex)  # interactions return的值: [atom_interaction, aa_interaction, key_aa]
        return interactions

    def metrics(self, true_interaction, pred_interactions):
        ious = [[iou_metrics(true_interaction[i], pred_interaction[i]) for i in range(3)] for pred_interaction in pred_interactions]
        screening_iou = [iou_virtual_screen_metrics(true_interaction[2], pred_interaction[2]) for pred_interaction in pred_interactions]
        atom_interact_iou, aa_interact_iou, key_aa_iou = list(zip(*ious))
        return [atom_interact_iou, aa_interact_iou, key_aa_iou, screening_iou]


def plt_iou(y, title):
    x = np.arange(0, len(y))
    y2 = np.array(y)
    plt.scatter(x, y2)
    plt.axhline(y=0.5, color='r', linestyle='-')
    plt.title(title)
    plt.show()


def plt_ious(datas):
    for y, title in datas:
        plt_iou(y, title)


# def clear_file():
#     os.system(f'rm -fr /dev/shm/__*')
#     os.system(f'rm -fr *.sdf')

def iou_compute(pair_data, parallel=1, plt_scatter=False, merge_complex='pymol'):
    MP = MultiProcessing2(merge_complex)
    pair_preds = [[pair_data[0][0], pair_pred[1]] if not pair_pred[0] else pair_pred for pair_pred in pair_data[1]]
    all_pair_data = [pair_data[0]] + pair_preds
    if parallel > 1:
        with Pool(parallel) as pool:
            all_interactions = pool.map(MP.get_interaction, all_pair_data)
    else:
        all_interactions = [MP.get_interaction(pair) for pair in all_pair_data]
    atom_interact_iou, aa_interact_iou, key_aa_iou, screening_iou = MP.metrics(all_interactions[0], all_interactions[1:])
    # 清除产生的中间文件
    # clear_file()  # 因为要分开到多个机器跑，所以有问题的

    if plt_scatter:
        plt_ious([(atom_interact_iou, '原子粒度IOU'), (aa_interact_iou, '残基粒度IOU'), (key_aa_iou, '关键残基IOU'),
                  (screening_iou, '关键残基复现比')])
    return atom_interact_iou, aa_interact_iou, key_aa_iou, screening_iou


if __name__ == '__main__':
    # pair_data = [(['example_data/1a30/1a30_p_true.pdb', 'example_data/1a30/1a30_l_true.sdf'], [None, 'example_data/1a30/1a30_l_pred.sdf'])]
    # iou_compute(pair_data, parallel=2, plt_scatter=False, merge_complex='pymol')
    print("Done!")
