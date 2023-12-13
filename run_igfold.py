#!/data/corp/hongning.zhang/software/igfold/0.4.0/.conda/bin/python3
# -*- coding: utf-8 -*-
from igfold import IgFoldRunner
from igfold.refine.pyrosetta_ref import init_pyrosetta
from pathlib import Path
from Bio import SeqIO
from multiprocessing import Pool
from typing import Dict

def run(
        sequences: Dict[str, str],
        output_file:Path,
        do_refine=True,
        do_renum=True,
        igfold: IgFoldRunner = None):

    init_pyrosetta()

    if igfold is None:
        igfold = IgFoldRunner()

    print(f"Predicting {output_file.stem}")
    return igfold.fold(
        str(output_file), # Output PDB file
        sequences=sequences, # Antibody sequences
        do_refine=do_refine, # Refine the antibody structure with PyRosetta
        do_renum=do_renum, # Renumber predicted antibody structure (Chothia)
    )


def predict_VHH(
          input_file: Path,
          output_dir: Path,
          do_refine=True,
          do_renum=True,
          n_processes: int = 1):
    seqs = SeqIO.parse(input_file, "fasta")
    output_dir.mkdir(parents=True, exist_ok=True)
    params = ((
         {"H": str(seq.seq)},
         output_dir / f"{seq.id}.pdb", 
         do_refine, do_renum) for seq in seqs)
    if n_processes > 1:
            with Pool(n_processes) as pool:
                pool.starmap(run, params)
    else:
        igfold = IgFoldRunner()
        for param in params:
            run(*param, igfold=igfold)


def predict_mAb(
          input_file: Path,
          output_dir: Path,
          do_refine=True,
          do_renum=True,
          n_processes: int = 1):
    seqs = SeqIO.parse(input_file, "fasta")
    output_dir.mkdir(parents=True, exist_ok=True)
    params = []
    for seq in seqs:
         output_file = output_dir / f"{seq.id}.pdb"
         heavy, light = str(seq.seq).split(":")
         params.append((
                {"H": heavy, "L": light},
                output_file, 
                do_refine, do_renum))

    if n_processes > 1:
            with Pool(n_processes) as pool:
                pool.starmap(run, params)
    else:
        igfold = IgFoldRunner()
        for param in params:
            run(*param, igfold=igfold)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(
         "input_file", type=Path, help="Input FASTA file")
    parser.add_argument(
         "output_dir", type=Path, help="Output directory")
    parser.add_argument(
         "--n_processes", '-n', type=int, default=1, help="Number of processes")
    parser.add_argument(
         "--do_refine", 
         action="store_true", 
         help="Refine the antibody structure with PyRosetta")
    parser.add_argument(
         "--do_renum", 
         action="store_true", 
         help="Renumber predicted antibody structure (Chothia)")
    parser.add_argument(
            "--format",
            choices=["VHH", "mAb"],
            default="VHH",
            help="Input antibody format")
    args = parser.parse_args()

    if args.format == "VHH":
        predict_func = predict_VHH
    elif args.format == "mAb":
        predict_func = predict_mAb
    else:
        raise ValueError(f"Unknown format {args.format}")

    predict_func(
        args.input_file,
        args.output_dir,
        args.do_refine,
        args.do_renum,
        args.n_processes)
