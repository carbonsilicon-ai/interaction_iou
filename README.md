## 1. 作用力评估指标定义 ###
### 1) 原子粒度IOU ### 
    蛋白配体复合物中，考虑蛋白原子索引，配体原子索引，相互作用力类型，得到原子粒度IOU。（protein_atom, ligand_atom, interaction_type）
### 2) 残基粒度IOU ### 
    复合物中，考虑蛋白残基，相互作用力类型，不考虑配体原子索引(对接多个不同分子导致无法考虑)，得到残基粒度IOU。（protein_res, interaction_type）
### 3) 关键残基IOU ###
    复合物中，考虑发生相互作用的蛋白残基，得到关键残基IOU。(protein_res)
### 4) 关键残基复现比 ###
    类似关键残基IOU,不同点在于交并比中的并只考虑晶体结构，以后不一定会用。只是先做了定义。
## 2. 使用方法及参数设置 ##
① 只计算iou,参考iou_test.py   
② 虚筛,参考iou_screening.py   
③ 虚筛后指标计算,参考forward_screening_power2.py(来源于沈超rtmscore)   

    atom_interact_iou, aa_interact_iou, key_aa_iou, screening_iou = iou_compute(pair_data, parallel=parallel, plt_scatter=plt_scatter)
    
    atom_interact_iou：以coreset数据集(285个数据)为例，所有数据的原子粒度IOU组成的元组。 
    aa_interact_iou：以coreset数据集(285个数据)为例，所有数据的残基粒度IOU组成的元组。
    key_aa_iou：以coreset数据集(285个数据)为例，所有数据的关键残基IOU组成的元组。
    screening_iou：以coreset数据集(285个数据)为例，所有关键残基复现比组成的元组。

### pair_data ###
准备好pair_data即可，tuple类型，第一个元素为晶体的蛋白配体文件对，第二个元素是所有预测的蛋白配体文件对组成的列表。  
一个tuple保证是同一个蛋白靶点下所有数据，多个蛋白靶点则for循环解决。单个蛋白靶点，具体格式如下。

    ([晶体蛋白1, 晶体配体1], [[晶体蛋白1, 预测配体1], [晶体蛋白1, 预测配体1'], [晶体蛋白1, 预测配体2], ...])
     
    
其中，由于"晶体蛋白1"出现多次，所以在预测的结果文件中的位置，可以为None(如果柔性对接，此处不能为None)。具体示例如下。

      (['example_data/true/1e66/1e66_protein.pdb', 'example_data/true/1e66/1e66_ligand.sdf'], [[None, 'example_data/pred/1e66/1e66_1a30/1a30_ligand_4.sdf'], 
                                                                                              [None, 'example_data/pred/1e66/1e66_1a30/1a30_ligand_5.sdf'], 
                                                                                              [None, 'example_data/pred/1e66/1e66_1bcu/1bcu_ligand_32.sdf'],
                                                                                              [None, 'example_data/pred/1e66/1e66_1bcu/1bcu_ligand_95.sdf']]) 
### parallel ###
    并行进程数，为1则不并行，默认为1。在腾讯服务器测试结果显示，parallel=10最佳。
### plt_scatter ###
    是否画散点图，默认不画
## 3. 格式问题 ##
    蛋白为pdb文件。
    晶体配体和预测配体格式一致，推荐都为sdf(每个sdf文件中只有一个构象。其余格式没测试过，可能会存在问题)。
    是否去氢?都可以。因为在pymol里会自动去氢，且plip官方复合物pdb文件里面配体会加极性氢(极性氢参与某些相互作用力的形成)。
### 错误示范 ###
    晶体配体格式为sdf，预测配体格式为pdb，虽然过程不报错，但是最终结果不对(原因在于pymol生成的复合物文件原子顺序不一致，所以保持晶体配体和预测配体均为sdf格式)。
### NOTE
    1. 需安装plip和pymol。
    2. 晶体配体和预测配体所在文件中的原子顺序需保持一致，否则计算错误。
    3. 相同数据多次计算结果有些许差距，plip软件本身算法导致，影响不大。
    4. 柔性对接(即蛋白在动)情况下，pair_data中不能为None。
    5. merge_complex = 'pymol',此处不可选，目前不支持prody。
       prody和pymol两种合成复合物的方式，实测结果指标差不多，默认pymol，调研大家用这个比较多。其中，prody这块使用了rdkit,可能会因为分子化合价等问题导致报错。目前没遇到用prody的必要，所以暂时没加入prody方式合成复合物的代码。
    


