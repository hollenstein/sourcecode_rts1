import os
import max_quant as MQ


# Process the MaxQuant results of the stress response and stress response Hog1
#   inhibition experiments
if __name__ == '__main__':
    root = '../'
    table_dir = os.path.join(root, 'tables')
    fasta_path = os.path.join(root, 'yeast_cont_20140324.fasta')
    mq_txt_folder = 'MQ_txt_SR_Romanov2017'
    evidence_path = os.path.join(
        root, 'MaxQuant', mq_txt_folder, 'evidence.txt'
    )
    experiments_toreverse = [
        '01_nacl_5min', '02_nacl_5min', '03_nacl_5min', '04_nacl_5min',
        '05_nacl_5min', '06_nacl_5min', '01_hog1as_5min', '02_hog1as_5min'
    ]
    setups = ['SR', 'SR-hog1as']
    
    from collections import defaultdict as ddict
    exp_to_setup = ddict(lambda: 'exclude')
    for exp in ['01_hog1as_5min', '02_hog1as_5min']:
        exp_to_setup[exp] = 'SR-hog1as'
    for exp in ['01_nacl_5min', '02_nacl_5min', '03_nacl_5min',
                '04_nacl_5min', '05_nacl_5min', '06_nacl_5min']:
        exp_to_setup[exp] = 'SR'

    # Import Evidence
    evidence, normalization_data = MQ.process_evidence(evidence_path, fasta_path)
    evidence['Setup'] = [exp_to_setup[exp] for exp in evidence['Experiment']]
    evidence = evidence[(evidence['Setup'] != 'exclude')]
    MQ.reverse_ratios(evidence, experiments_toreverse)

    # Generate output files
    for setup in setups:
        ev_setup = evidence[(evidence['Setup'] == setup)]

        evidence_outpath = os.path.join(table_dir, 'evidence_' + setup + '.tsv')
        phosphosite_outpath = os.path.join(table_dir, 'phospho_' + setup + '.tsv')
        protein_outpath = os.path.join(table_dir, 'protein_' + setup + '.tsv')

        MQ.write_evidence_table(ev_setup, evidence_outpath)
        MQ.write_phosphosite_table(ev_setup, phosphosite_outpath, prob_cutoff=0.7)
        MQ.write_protein_table(ev_setup, protein_outpath)
