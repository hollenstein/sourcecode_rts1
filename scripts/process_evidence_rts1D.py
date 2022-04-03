import os
import pandas as pd
import max_quant as MQ


# Process the MaxQuant results of the rts1 delta experiment
if __name__ == '__main__':
    root = '../'
    table_dir = os.path.join(root, 'tables')
    fasta_path = os.path.join(root, 'yeast_cont_20140324.fasta')
    mq_txt_folder1 = 'MQ_txt_rts1D_1'
    mq_txt_folder2 = 'MQ_txt_rts1D_2'
    evidence_path1 = os.path.join(
        root, 'MaxQuant', mq_txt_folder1, 'evidence.txt'
    )
    evidence_path2 = os.path.join(
        root, 'MaxQuant', mq_txt_folder2, 'evidence.txt'
    )
    experiments_toreverse = ['04_rts1', 'F_rts1', 'I_rts1', '08_rts1_1']
    setup = 'rts1D'

    # Import Evidence
    evidence1, normalization_data = MQ.process_evidence(evidence_path1, fasta_path)
    evidence2, normalization_data = MQ.process_evidence(evidence_path2, fasta_path)
    evidence = pd.concat([evidence1, evidence2])
    evidence['Setup'] = setup
    MQ.reverse_ratios(evidence, experiments_toreverse)

    # Generate output files
    evidence_outpath = os.path.join(table_dir, 'evidence_' + setup + '.tsv')
    phosphosite_outpath = os.path.join(table_dir, 'phospho_' + setup + '.tsv')
    protein_outpath = os.path.join(table_dir, 'protein_' + setup + '.tsv')

    MQ.write_evidence_table(evidence, evidence_outpath)
    MQ.write_phosphosite_table(evidence, phosphosite_outpath, prob_cutoff=0.7)
    MQ.write_protein_table(evidence, protein_outpath)
