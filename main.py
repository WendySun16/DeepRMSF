import numpy as np
from moleculekit.molecule import Molecule
import os
import torch
import torch.nn as nn
import subprocess
import argparse
import sys
sys.path.append('..')
import mrcfile
import string
import shutil
from rna_util import get_dssr, get_file, rna_pdb2sse, rna_get_smi_map, select_a, rna_get_dir
from rna_to_input import rna_save_ana_map, rna_get_box_list, rna_get_sp_data
from model_util import rna_concat_data
from rna_five_fold import test_model_for_rmsf
from rmsf_model import rmsf_model, init_weights


parser = argparse.ArgumentParser(description='DeepRMSF')
parser.add_argument('--ori_dir', type=str, default="/DeepRMSF/test", help='folder for saving PDBs of RNAs')
parser.add_argument('--data_dir', type=str, default="/data/map_data", help='folder where you want to save the simulated maps')
parser.add_argument('--box_file', type=str, default="/data/map_data/datafile", help='folder where you want to save the box_files')
parser.add_argument('--log_dir', type=str, default="/DeepRMSF/test", help='folder where you want to save predicted data')

args = parser.parse_args()


def main():
    # step one: use x3dna-dssr
    ori_dir = args.ori_dir
    for root, dirs, files in os.walk(ori_dir):
        for dir in dirs:
            get_dssr(dir, ori_dir)
            org_pdb = f"{ori_dir}/{dir}/{dir}.pdb"
            # step two：determine the secondary structure of each residue and generate a file
            dssr_pairs = get_file(dir, ori_dir)
            rna_pdb2sse(org_pdb, dssr_pairs)
            # step three：simulated map
            fin_pdb = f"{ori_dir}/{dir}/{dir}_a.pdb"
            select_a(org_pdb, fin_pdb)
            out_file = f"{ori_dir}/{dir}/{dir}_out.mrc"
            output = rna_get_smi_map(fin_pdb, out_file)
            print(output)

    data_dir = args.data_dir
    pdbid_list = rna_get_dir(ori_dir)
    for pdbid in pdbid_list:
        if not os.path.exists(f"{data_dir}/{pdbid}"):
            os.mkdir(f"{data_dir}/{pdbid}")
    rna_save_ana_map(pdbid_list, data_dir, ori_dir, pre=True)

    box_file = args.box_file
    for pdbid in pdbid_list:
        out_file = f"{box_file}/{pdbid}.pth"
        rna_get_box_list(data_dir, [pdbid], out_file, sim=True, pre=True, res=4)

    test_list = []
    for pdbid in pdbid_list:
        test_list.append(pdbid + ".pth")
    test_data = rna_concat_data(box_file, test_list, dsize=(40,40,40,40,10), train=False, only_rmsf=False, only_sse=False)
    torch.save(test_data, f"{ori_dir}/test_data.pth")
    
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ['CUDA_VISIBLE_DEVICES'] ='0'
    log_dir = args.log_dir
    record = f"{log_dir}/record.txt"
    # initialize the model, and then call the trained parameters
    model = rmsf_model(in_channels = 424)
    init_weights(model)
    model_file = "/data/model/5/model.pth"
    state_dict = torch.load(model_file)
    model_state_dict = model.state_dict()
    for key in state_dict:
        if key[7:] in model_state_dict:
            model_state_dict[key[7:]] = state_dict[key]
    model.load_state_dict(model_state_dict)
    with open(record,'w') as f:
        _, time_list = test_model_for_rmsf(model, 5, log_dir, ori_dir, ori_dir, f, sse=True, pre=True)

    # output dot-bracket
    path = f"{ori_dir}/{pdbid}"
    os.chdir(path)
    shutil.copy("dssr-2ndstrs.dbn", log_dir)


if __name__ == '__main__':
    main()