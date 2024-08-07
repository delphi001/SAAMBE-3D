#!/usr/bin/env python
# coding: utf-8


import os
import sys
import warnings
import argparse
import numpy as np
import prody as pdy
import xgboost as xgb

from utils.protseqfeature import *

__version__ = "1.0"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Let's suppress the log messages from ProDy, to keep stdout clear.
pdy.confProDy(verbosity='none')
# Let's suppress non-critical XGBoost warnings arising because model is developed 
# using an older verion 0.87.
xgb.set_config(verbosity=0)


def mutation_pdb(pdbid):
    import urllib.request, urllib.parse, urllib.error

    url = "https://files.rcsb.org/view/" + pdbid + ".pdb"
    urllib.request.urlretrieve(url, pdbid)


def mutation_sequence(pdbid, resid, chain):
    label_index = 1
    resid_label = []
    resid_label_aa = []
    mutation_coordinate = []
    resid_label_aa = [0 for _ in range(11)]  # how many sequence to use
    print(">>mutation_sequence >> ", pdbid, resid, chain)
    for i in range(len(resid_label_aa)):
        resid_int = "".join([ch for ch in resid.strip() if ch in "-0123456789"])
        resid_label.append(int(resid_int) - (len(resid_label_aa) - 1) // 2 + i)
        # print("resid_label>>", resid_label)
        for line in open(pdbid):
            pdbstr = line.strip()
            if (
                pdbstr[0:4] == "ATOM"
                and pdbstr[21:22] == chain
                and pdbstr[13:15] == "CA"
            ):
                if pdbstr[22:27].strip() == str(resid_label[i]) or (
                    pdbstr[22:27].strip() == str(resid) and label_index == 1
                ):
                    if pdbstr[22:27].strip() == str(resid):
                        mutation_coordinate = [
                            float(pdbstr[29:38].strip()),
                            float(pdbstr[38:46].strip()),
                            float(pdbstr[46:55].strip()),
                        ]
                        label_index = 0
                    resid_label_aa[i] = pdbstr[17:20]
                    # print("pdbstr>>",pdbstr[13:15], pdbstr[17:20], pdbstr[22:27].strip(), str(resid_label[i]))
                    break
    if len(mutation_coordinate) != 3:
        error_index = 1
    else:
        error_index = 0
    return resid_label_aa, mutation_coordinate, error_index


def mutation_distance(pdbid, chain, mutation_coordinate):
    resid_label_dis_aa = []
    resid_label_aa = []
    resid_label_distance = []
    for line in open(pdbid):
        pdbstr = line.strip()
        if pdbstr[0:4] == "ATOM" and pdbstr[13:15] == "CA" and pdbstr[21:22] != chain:
            mutation_coordinate1 = [
                float(pdbstr[29:38].strip()),
                float(pdbstr[38:46].strip()),
                float(pdbstr[46:55].strip()),
            ]
            resid_label_aa.append(pdbstr[17:20])
            resid_label_distance.append(
                np.sqrt(
                    np.square(mutation_coordinate[0] - mutation_coordinate1[0])
                    + np.square(mutation_coordinate[1] - mutation_coordinate1[1])
                    + np.square(mutation_coordinate[2] - mutation_coordinate1[2])
                )
            )
    b = list(zip(resid_label_distance, list(range(len(resid_label_distance)))))
    b.sort(key=lambda x: x[0])
    c = [x[1] for x in b]
    sequ_num = 10  # 10 sequence in total,use last 7sequence
    if len(c) >= sequ_num:
        for j in range(0, sequ_num, 1):
            if b[j][0] <= 10:
                resid_label_dis_aa.append(resid_label_aa[c[j]])
            else:
                resid_label_dis_aa.append(0)
    else:
        for j in range(0, len(c), 1):
            if b[j][0] <= 10:
                resid_label_dis_aa.append(resid_label_aa[c[j]])
            else:
                resid_label_dis_aa.append(0)
        while len(resid_label_dis_aa) < sequ_num:
            resid_label_dis_aa.append(0)
    return resid_label_dis_aa


def mutation_pdb_information(pdbid):
    resolution, r_value, temperature, ph = 0, 0, 0, 0
    for line in open(pdbid):
        pdbstr = line.strip()
        if pdbstr[0:22] == "REMARK   2 RESOLUTION.":
            try:
                resolution = float(pdbstr[26:30].strip())
            except ValueError:
                resolution = 0
        if pdbstr[0:45] == "REMARK   3   R VALUE            (WORKING SET)":
            try:
                r_value = float(pdbstr[49:54].strip())
            except ValueError:
                r_value = 0
        if pdbstr[0:23] == "REMARK 200  TEMPERATURE":
            try:
                temperature = float(pdbstr[45:48].strip())
            except ValueError:
                temperature = 0
        if pdbstr[0:14] == "REMARK 200  PH":
            ph = "".join([ch for ch in pdbstr[45:48].strip() if ch in "0123456789."])
            try:
                ph = float(ph)
            except ValueError:
                ph = 0
            break
    return resolution, r_value, temperature, ph


def mutation_sequence2(struct, resid, chain, window=(-5, 5)):
    ch = chain
    resi = int("".join([c for c in resid if c in "-0123456789"]))
    wt_sel_str = f"chain {ch} and resid {resi} and name CA"
    wt_sel = struct.select(wt_sel_str)
    resid_label_aa = [0 for _ in range(window[0], window[1] + 1)]
    if not (wt_sel is None) and len(wt_sel.getResnames()) == 1:
        mt_coords = wt_sel.getCoords()[0]
        chain_ca = struct.select(f"chain {ch} and name CA")
        chain_resi = [
            (ch, ri, rn)
            for ri, rn in zip(chain_ca.getResnums(), chain_ca.getResnames())
        ]
        ri_list = [(i, e) for i, e in enumerate(chain_resi) if e[1] == resi]
        if len(ri_list) == 1:
            rindx = ri_list[0][0]
            lbl_shift = 0
            rindx_from, rindx_to = rindx + window[0], rindx + window[1] + 1
            if rindx_to > len(chain_resi):
                rindx_to = len(chain_resi)
            if window[1] - rindx > 0:
                rindx_from = 0
                lbl_shift = window[1] - rindx
            for i, e in enumerate(chain_resi[rindx_from:rindx_to]):
                resid_label_aa[lbl_shift + i] = e[2]
                
    return resid_label_aa, mt_coords


def validate_mt_info(struct, chain, resid, wt, mt, window=(-5, 5)):
    warns = []
    errors = []
    ch = chain
    wt_sel_str = f"chain {ch} and resid {resid} and name CA"
    wt_sel = struct.select(wt_sel_str)
    if not (wt_sel is None) and len(wt_sel.getResnames()) == 1:
        resn = aa3to1.get(wt_sel.getResnames()[0].upper(), "X")
        # neigh_resn = [aa3to1.get(aa3.upper(), 'X') for aa3 in neigh_sel.getResnames()]
        # print(neigh_resn, resn)
        if resn.upper() != wt.upper():
            errors.append(
                f"ERROR>> Claimed wild-type residue '{wt}' in 'chain {ch} and resid {resid}' found '{wt_sel.getResnames()[0]}'."
            )
        resi = int("".join([c for c in resid if c in "-0123456789"]))
        chain_ca = struct.select(f"chain {ch} and name CA")
        chain_resi = [
            (ch, ri, rn)
            for ri, rn in zip(chain_ca.getResnums(), chain_ca.getResnames())
        ]
        ri_list = [(i, e) for i, e in enumerate(chain_resi) if e[1] == resi]
        if len(ri_list) == 1:
            rindx = ri_list[0][0]
            lbl_shift = 0
            rindx_from, rindx_to = rindx + window[0], rindx + window[1] + 1
            if rindx_to > len(chain_resi):
                rindx_to = len(chain_resi)
                warns.append(
                    f"WARNING>> Missing some sequence neighbors 'chain {ch} and resid {resid}' of +5 positions"
                )
            if window[1] - rindx > 0:
                rindx_from = 0
                lbl_shift = window[1] - rindx
                warns.append(
                    f"WARNING>> Missing some sequence neighbors 'chain {ch} and resid {resid}' of -5 positions"
                )
        neigh_sel = struct.select(
            f"(within 15 of ({wt_sel_str})) and (name CA and not chain {ch})"
        )
        if not neigh_sel is None:
            for ci, ri, aa3 in zip(
                neigh_sel.getChids(), neigh_sel.getResnums(), neigh_sel.getResnames()
            ):
                if not aa3 in aa3to1:
                    warns.append(
                        f"WARNING>> In vicinity of 'chain {ch} and resid {resid}' found non-standard residue 'chain {ci} and resid {ri} and resn {aa3}'."
                    )
    elif wt_sel is None:
        errors.append(
            f"ERROR>> Claimed wild-type residue '{wt}' in 'chain {ch} and resid {resid}' not found."
        )
    return warns, errors


def check_pdb(pdbfile, aa3to1, mt_info=None, mt_file=None):
    warns_list = []
    errors_list = []
    struct = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        struct = pdy.parsePDB(pdbfile, header=False)

    if len(struct.getCoordsets()) > 1:
        warns_list.append(
            "WARNING>> Multiple models detected in the input pdb file, only first model coordinated will be used."
        )
    if len(set(struct.getAltlocs())) > 1:
        warns_list.append(
            "WARNING>> Alternate locations detected for some residues in the input pdb file, "
            "only highest occupancy coordinates for alternate locations of atoms will be utilized."
        )
    if (not mt_info is None) and (len(mt_info) == 4):
        warns, errors = validate_mt_info(
            struct, mt_info[0], mt_info[1], mt_info[2], mt_info[3]
        )
        if len(warns) > 0:
            warns_list.extend(warns)
        if len(errors) > 0:
            errors_list.extend(errors)
    elif mt_file:
        for ln in open(mt_file):
            if ln.strip().startswith("#") or len(ln.strip()) == 0:
                continue
            ws = ln.strip().split()
            if len(ws) == 4:
                warns, errors = validate_mt_info(struct, ws[0], ws[1], ws[2], ws[3])
                if len(warns) > 0:
                    warns_list.extend(warns)
                if len(errors) > 0:
                    errors_list.extend(errors)
            else:
                warns_list.append((
                    f"ERROR>> Mutation list file format error, expected "
                    f"four words, chain resid wt mt but found: '{ln.strip()}'."
                ))
    if len(warns_list) > 0:
        print("\n".join(warns_list), file=sys.stderr)
    if len(errors_list) > 0:
        print("\n".join(errors_list), file=sys.stderr)
    return struct, warns_list, errors_list


def file_loop(
    pdb_id, struct, mutation_chain, mutation_resid, wild_aa, mutation_aa, model_type, verbose
):
    if str(model_type) == "1":
        if wild_aa == mutation_aa:
            return "0 Neutral"
    else:
        if wild_aa == mutation_aa:
            return "Neutral"
    if len(wild_aa) == 3:
        wild_aa = aa3to1.get(wild_aa, 'X')
    if len(mutation_aa) == 3:
        mutation_aa = aa3to1.get(mutation_aa, 'X')
    label = []
    label.append(net_volume2(aa1_map, wild_aa, mutation_aa))
    label.append(net_hydrophobicity2(aa1_map, wild_aa, mutation_aa))
    label.append(net_flexibility2(aa1_map, wild_aa, mutation_aa))
    label.append(mutation_hydrophobicity2(aa1_map, wild_aa, mutation_aa))
    label.append(mutation_polarity2(aa1_map, wild_aa, mutation_aa))
    label.append(mutation_type2(aa1_mttype_label, wild_aa, mutation_aa))
    label.append(mutation_size2(aa1_map, wild_aa, mutation_aa))
    label.append(mutation_hbonds2(aa1_map, wild_aa, mutation_aa))
    label.append(mutation_chemical2(aa1_map, wild_aa, mutation_aa))
    if os.path.isfile(pdb_id) == 0:
        return "Please check the PDB file"
    (resid_label_aa, mutation_coordinate) = mutation_sequence2(
        struct, str(mutation_resid), mutation_chain
    )
    
    label_aa_distance = mutation_distance(pdb_id, mutation_chain, mutation_coordinate)
    for i in resid_label_aa:
        label.append(aa3_label.get(i, 0))
    for i in label_aa_distance:
        label.append(aa3_label.get(i, 0))
    (reso, r_value, temp, ph) = mutation_pdb_information(pdb_id)
    label.append(reso)
    label.append(temp)
    label.append(ph)
    if verbose:
        print(
            "Mutation features",
            pdb_id,
            mutation_chain,
            mutation_resid,
            wild_aa,
            mutation_aa,
            model_type,
            label,
        )

    return pred_feature(label, model_type)


def pred_feature(label, model_type):
    if str(model_type) == "1":
        model = xgb.Booster(model_file=f"{__location__}/regression_v01.model")
        # model.save_model(f"{__location__}/regression_v01.model")
    else:
        model = xgb.Booster(model_file=f"{__location__}/classification_v01.model")
        # model.save_model(f"{__location__}/classification_v01.model")
    x = np.array(label, dtype=object)
    x = x.reshape((1, len(label)))
    x = xgb.DMatrix(x)
    y_pred = model.predict(x)
    if str(model_type) == "1":
        if y_pred[0] > 0:
            return "%.2f" % y_pred[0] + " Destabilizing"
        else:
            return "%.2f" % y_pred[0] + " Stabilizing"
    else:
        if y_pred[0] > 0.5:
            return "Disruptive"
        else:
            return "Nondisruptive"


def output_format(model):
    if str(model) == "1":
        return "Structure_file Chain Position Wild Mutant ddG(kcal/mol) Type"
    else:
        return "Structure_file Chain Position Wild Mutant Type"


def is_file(filepath):
    val = os.path.isfile(filepath)
    if not val:
        raise argparse.ArgumentTypeError("`%s` is not a valid filepath" % filepath)
    return filepath


def argument_parser():
    parser = argparse.ArgumentParser(
        prog="saambe-3d.py",
        description=(f"SAAMBE-3D_v{__version__} : "
                     "Predict the free energy change of binding due to "
                     "point mutation for protein-protein binding"),
    )
    parser.add_argument(
        "-i",
        "--complex-pdb",
        help="PDB file of the wild-type protein-protein complexz",
        type=is_file,
        required=True,
    )
    parser.add_argument(
        "-d",
        "--model-type",
        type=int,
        help="model type: 1. prediction of ddG, 0. prediction of stabilizing/destabilizing. (default: 1)",
        choices={1, 0},
        default=1,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output file name for the predictions (default: output.out)",
        default="output.out",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="print verbose messages",
        action="store_true",
    )
    group1 = parser.add_argument_group(
        "single point-mutation", "single point-mutation information"
    )
    group1.add_argument(
        "-c", "--chain", help="chain of the mutation in the input PDB file"
    )
    group1.add_argument(
        "-r", "--resid", help="resid of mutation position",
    )
    group1.add_argument(
        "-w", "--wild-type", help="wild-type amino acid's 1 letter code",
    )
    group1.add_argument(
        "-m", "--mutant", help="mutant amino acid's 1 letter code",
    )
    group2 = parser.add_argument_group(
        "multiple point-mutations", "multiple point-mutations listing"
    )
    group2.add_argument(
        "-f",
        "--mutation-list-file",
        help="".join(
            [
                "Mutation list file. The file should have one mutation per line ",
                "where every line following four information seperated by space. ",
                "chain resid wild-type(1-letter-code) mutant(1-letter-code)",
            ]
        ),
    )
    args = parser.parse_args()
    return args


def main():
    args = argument_parser()
    struct = None
    if args.mutation_list_file:
        f = open(args.output, "w")
        print(output_format(args.model_type), file=f)
        struct, warns_list, errors_list = check_pdb(
            args.complex_pdb, aa3to1, mt_file=args.mutation_list_file
        )
        if len(errors_list) > 0:
            sys.exit()
        for line in open(args.mutation_list_file):
            if line.strip().startswith("#"):
                continue
            info = line.strip().split(" ")
            pred = file_loop(
                args.complex_pdb,
                struct,
                info[0],
                info[1],
                str(info[2]),
                info[3],
                args.model_type,
                args.verbose
            )
            print(
                args.complex_pdb, info[0], info[1], str(info[2]), info[3], pred, file=f
            )
        f.close()
        print(f"Check outputs in file: '{args.output}'")
    else:
        f = open(args.output, "w")
        struct, warns_list, errors_list = check_pdb(
            args.complex_pdb,
            aa3to1,
            (args.chain, args.resid, args.wild_type, args.mutant),
        )
        if len(errors_list) > 0:
            sys.exit()
        pred = file_loop(
            args.complex_pdb,
            struct,
            args.chain,
            args.resid,
            args.wild_type,
            args.mutant,
            args.model_type,
            args.verbose
        )
        print(pred, file=f)
        f.close()
        print(f"Check outputs in file: '{args.output}'")


if __name__ == "__main__":
    main()
