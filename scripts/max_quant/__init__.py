from collections import OrderedDict as odict
import numpy as np
import pandas as pd
import maspy.peptidemethods
import maspy.proteindb


def write_evidence_table(ev, path):
    """ Write an updated evidence table file. """
    columns = [
        'Sequence', 'Modified sequence', 'Sequence start', 'Sequence end',
        'Protein Systematic Name', 'Protein Standard Name', 'Phospho (STY)',
        'Phospho positions', 'Phospho sites', 'Protein phospho sites',
        'Ratio Log2 normalized', 'Ratio Log2', 'Combined Phospho Probability',
        'Valid Quantification', 'Raw file', 'Experiment', 'Setup',
        'Reverse', 'Potential contaminant', 'S/T-P motif',
    ]
    ev.ix[:, columns].to_csv(path, sep='\t', index=False)


def write_phosphosite_table(ev, path, prob_cutoff=0):
    """ Write a phosphosite table. """
    m = np.all([
        ev['Valid Quantification'],
        ev['Phospho (STY)'] > 0,
        ev['Combined Phospho Probability'] >= prob_cutoff,
    ], axis=0)

    headers = [
        'Protein Standard Name', 'Protein Systematic Name',
        'Phospho (STY)', 'Ratio Log2 normalized',
        'Phospho positions', 'Phospho sites', 'S/T-P motif', 'Setup'
    ]
    aggregation = odict([(h, 'first') for h in headers])
    aggregation['Ratio Log2 normalized'] = 'mean'
    groupby = ['Experiment', 'Protein phospho sites']
    phospho_table = ev[m].groupby(groupby, as_index=False).agg(aggregation)
    phospho_table.to_csv(path, sep='\t', index=False)


def write_protein_table(ev, path):
    """ Write a protein table. """
    m = ev['Valid Quantification'] & (ev['Phospho (STY)'] == 0)

    headers = ['Protein Systematic Name', 'Ratio Log2 normalized', 'Setup']
    aggregation = odict([(h, 'first') for h in headers])
    aggregation['Ratio Log2 normalized'] = 'mean'
    groupby = ['Experiment', 'Protein Standard Name']
    protein_table = ev[m].groupby(groupby, as_index=False).agg(aggregation)
    protein_table.to_csv(path, sep='\t', index=False)


def process_evidence(evidence, fasta):
    """ Read and process a MaxQuant evidence file. """
    ev = pd.read_csv(evidence, sep='\t')
    _adjust_mq_version_headers(ev)
    proteindb = import_mq_proteindb(fasta)

    _add_contamination(ev, '[cont]')
    _add_decoy(ev, 'REV__')
    _add_proline_count(ev)
    _add_protein_names(ev, proteindb)
    _add_peptide_positions(ev, proteindb)
    _change_mod_sequence(ev)
    _add_protein_phospho_annotation(ev, proteindb)
    _add_stp_motif(ev, proteindb)
    _add_combined_probability(ev)

    _add_valid_quant(ev)
    _add_log_ratio(ev)
    ev, normalization_data = _normalize_ratios(ev)
    return ev, normalization_data


def reverse_ratios(ev, experiments):
    """ Invert log transformed ratios of the specified experiments. """
    for exp in experiments:
        m = (ev['Experiment'] == exp) & ev['Valid Quantification']
        ev.loc[m, 'Ratio Log2 normalized'] *= -1


def normalize_protein_abundance(ev, protein):
    """ Normalize peptide log ratios of protein by the protein abundance. 

    Peptides are normalized by subtracting the median log ratio of all
    non-phosphorylated peptides.
    """
    for exp in np.unique(ev['Experiment']):
        target_mask = np.all([
            ev['Valid Quantification'], ev['Experiment'] == exp,
            ev['Protein Standard Name'] == protein,
        ], axis=0)
        norm_mask = (ev['Phospho (STY)'] == 0) & target_mask
        median = np.median(ev['Ratio Log2 normalized'][norm_mask])
        ev.loc[target_mask, 'Ratio Log2 normalized'] -= median


def import_mq_proteindb(fastapath):
    proteindb = maspy.proteindb.importProteinDatabase(
        fastapath, minLength=4, maxLength=50, missedCleavage=10,
        ignoreIsoleucine=True, headerParser=maspy.proteindb.fastaParseSgd
    )
    # In contrary to MasPy, MaxQuant ignores Leucine and replaces "L" with "I"
    for peptide in list(proteindb.peptides):
        new_peptide = peptide.replace('L', 'I')
        if new_peptide not in proteindb.peptides:
            proteindb.peptides[new_peptide] = proteindb.peptides[peptide]
    # Add protein Standard Name EPO1 for YMR124W
    proteindb['YMR124W'].name = 'EPO1'
    return proteindb


def _add_protein_phospho_annotation(ev, proteindb):
    m = (ev['Phospho (STY)'] != 0) & (ev['Reverse'] != '+')
    all_phospho_positions = []
    all_phospho_sites = []
    all_protein_phospho_sites = []
    for row_id, row in ev[m].iterrows():
        protein = row['Leading Razor Protein']
        protein_std = proteindb.proteins[protein].name
        start, end = proteindb.peptides[row['Sequence']].proteinPositions[protein]
        mod_positions = maspy.peptidemethods.returnModPositions(
            row['Modified sequence'], indexStart=start
        )
        phospho_site_list = list()
        for phospho_pos in mod_positions['ph']:
            aa = proteindb.proteins[protein].sequence[phospho_pos - 1]
            site_string = str(phospho_pos) + '(' + aa + ')'
            phospho_site_list.append(site_string)

        phospho_positions = ','.join([str(i) for i in mod_positions['ph']])
        all_phospho_positions.append(phospho_positions)
        phospho_sites = ' / '.join(phospho_site_list)
        all_phospho_sites.append(phospho_sites)
        protein_phospho_sites = ' - '.join([protein_std, phospho_sites])
        all_protein_phospho_sites.append(protein_phospho_sites)

    ev['Phospho positions'] = ''
    ev.loc[m, 'Phospho positions'] = all_phospho_positions
    ev['Phospho sites'] = ''
    ev.loc[m, 'Phospho sites'] = all_phospho_sites
    ev['Protein phospho sites'] = ''
    ev.loc[m, 'Protein phospho sites'] = all_protein_phospho_sites


def _add_stp_motif(ev, proteindb):
    m = (ev['Phospho (STY)'] > 0) & (ev['Reverse'] != '+')
    stp_annotation = []
    for prot, sites in zip(
            ev[m]['Leading Razor Protein'],
            ev[m]['Phospho positions']
        ):
        annotation = []
        for site in map(int, sites.split(',')):
            motif = proteindb.proteins[prot].sequence[site - 1:site + 1]
            try:
                if motif in ['SP', 'TP']:
                    annotation.append(motif)
                else:
                    annotation.append('-')
            except IndexError:
                annotation.append('-')
        stp_annotation.append(' / '.join(annotation))
    ev['S/T-P motif'] = ''
    ev.loc[m, 'S/T-P motif'] = stp_annotation


def _change_mod_sequence(ev):
    func = lambda s: s.split('_')[1].replace('(', '[').replace(')', ']')
    ev['Modified sequence'] = ev['Modified sequence'].apply(func)


def _add_proline_count(ev):
    proline_count = [sequence.count('P') for sequence in ev['Sequence']]
    ev['P Count'] = proline_count


def _add_protein_names(ev, proteindb):
    m = (ev['Reverse'] != '+')
    protein_names = []
    for protein in ev['Leading Razor Protein'][m]:
        protein_names.append(proteindb[protein].name)

    ev['Protein Systematic Name'] = ev['Leading Razor Protein'].tolist()
    ev['Protein Standard Name'] = ev['Leading Razor Protein'].tolist()
    ev.loc[m, 'Protein Standard Name'] = protein_names


def _add_peptide_positions(ev, proteindb):
    m = (ev['Reverse'] != '+')
    start_positions = []
    end_positions = []
    for prot, seq in zip(ev['Leading Razor Protein'][m], ev['Sequence'][m]):
        start, end = proteindb.peptides[seq].proteinPositions[prot]
        start_positions.append(start)
        end_positions.append(end)
    ev['Sequence start'] = ''
    ev.loc[m, 'Sequence start'] = start_positions
    ev['Sequence end'] = ''
    ev.loc[m, 'Sequence end'] = end_positions


def _add_contamination(ev, contamination_tag):
    contamination = []
    for protein in ev['Leading Razor Protein']:
        if protein.find(contamination_tag) != -1:
            contamination.append('+')
        else:
            contamination.append('')
    ev['Potential contaminant'] = contamination


def _add_decoy(ev, decoy_tag):
    decoy = []
    for protein in ev['Leading Razor Protein']:
        if protein.find(decoy_tag) != -1:
            decoy.append('+')
        else:
            decoy.append('')
    ev['Reverse'] = decoy


def _add_valid_quant(ev):
    """ True if quantified and neiter a contaminant or a decoy protein. """
    valid_quantification = np.all([
        np.isfinite(ev['Ratio H/L']),
        ev['Potential contaminant'] != '+',
        ev['Reverse'] != '+',
    ], axis=0)
    ev['Valid Quantification'] = valid_quantification


def _add_combined_probability(ev):
    ev['Combined Phospho Probability'] = np.nan
    m = (ev['Phospho (STY)'] != 0) & (ev['Type'] != 'MULTI-MATCH')
    combined_probabilities = []
    for prob, num in zip(
            ev['Phospho (STY) Probabilities'][m],
            ev['Phospho (STY)'][m]
        ):
        combined_probability = _combined_localization_probability(prob, num)
        combined_probabilities.append(combined_probability)
    ev.loc[m, 'Combined Phospho Probability'] = combined_probabilities


def _combined_localization_probability(entry, num):
    probabilities = _extract_phospho_probabilities(entry)
    combined_probability = 1
    for probability, site in probabilities[:num]:
        combined_probability = combined_probability * probability
    return combined_probability


def _extract_phospho_probabilities(entry):
    """
    :params entry: entry of the field "Phospho (STY) Probabilities" from a
        maxquant evidence.txt file. For exaample: "AAADAIS(1)DIEIK"
    :returns: a sorted list (best probabilities first) of tuples, containing a
        probability value and the phospho position (starting from 1).
        For example: [(1.0, 7)]
    """
    probabilities = list()
    while entry.find('(') != -1:
        mod_position = int(entry.find('('))
        mod_probability = float(entry.split('(')[1].split(')')[0])
        probabilities.append((mod_probability, mod_position))
        entry = entry.split('(', 1)[0] + entry.split(')', 1)[1]
    return sorted(probabilities, reverse=True)


def _add_log_ratio(ev):
    m = ev['Valid Quantification']
    ev['Ratio Log2'] = np.nan
    ev.loc[m, 'Ratio Log2'] = np.log2(ev['Ratio H/L'][m])


def _normalize_ratios(ev):
    groups = []
    normalization_data = {}
    for exp, group in ev.groupby('Experiment'):
        norm_data = _calc_normalization(group)
        normalization_data[exp] = norm_data
        m = group['Valid Quantification']
        corr = (group['P Count'][m] * norm_data['p_shift']) + norm_data['shift']
        group['Ratio Log2 normalized'] = np.nan
        group['Ratio Log2 normalized'][m] = group['Ratio Log2'][m] - corr
        groups.append(group)
    return pd.concat(groups), normalization_data


def _calc_normalization(ev):
    p_masks = []
    for pcount in range(3):
        p_mask = np.all([
            ev['Valid Quantification'],
            ev['Phospho (STY)'] == 0,
            ev['P Count'] == pcount,
        ], axis=0)
        p_masks.append(p_mask)

    median_P0 = np.median(np.median(ev['Ratio Log2'][p_masks[0]]))
    median_P1 = np.median(np.median(ev['Ratio Log2'][p_masks[1]]))
    median_P2 = np.median(np.median(ev['Ratio Log2'][p_masks[2]]))

    weighted_p1 = (median_P1 - median_P0) / 1 * p_masks[1].sum()
    weighted_p2 = (median_P2 - median_P0) / 2 * p_masks[2].sum()
    total_num = p_masks[1].sum() + p_masks[2].sum()
    proline_shift = (weighted_p1 + weighted_p2) / total_num
    ratio_shift = median_P0

    norm_data = {'shift': ratio_shift, 'p_shift': proline_shift}
    return norm_data


def _adjust_mq_version_headers(ev):
    """ Unify evidence table headers of different MQ versions. """
    ev.rename(columns={
        'Leading razor protein': 'Leading Razor Protein',
        'Leading proteins': 'Leading Proteins'
    }, inplace=True)
