# Copyright (C) 2016 The President and Fellows of Harvard College
#
# This file may require additional software or modules to be installed to run
# run properly.
#
# The Genome Recoder software is available under an internal non-commercial
# research and academic use license.  Questions about this software or the
# licensing thereof can be addressed to  Office of Technology Development,
# Harvard University, email: otd@harvard.edu.
#
# @author Gleb Kuznetsov (kuznetsov@g.harvard.edu)

"""
Analysis scripts to identify patterns before recoding.
"""

import copy
import csv
import math
import os
import re

from Bio.Seq import reverse_complement
import pandas as pd

from getk.biopython_util import get_feature_by_gene_name
from getk.biopython_util import get_feature_by_id
from getk.biopython_util import get_feature_gene_name
from getk.biopython_util import get_feature_qualifier_value
from getk.biopython_util import get_genome_record
from getk.scorer import SSScorerOutsideFeature
from getk.scorer import PreserveRBSScorer
from getk.tool_wrappers.rbs_calc_util import calc_internal_rbs_expression
from getk.tool_wrappers.rbs_calc_util import calc_rbs_score_for_feature
from getk.tool_wrappers.rbs_calc_util import UNKNOWN_RBS_CALC_RESULT
from replace_agr import _build_refactor_context
from replace_agr import _check_recoding
from replace_agr import _get_original_genome_record
from replace_agr import MIX_NON_SYNONYMOUS_RBS_SS_0_65
from replace_agr import MIX_RBS_SS_0_2
from replace_agr import MIX_RBS_SS_0_5
from replace_agr import MIX_RBS_SS_0_65
from replace_agr import MIX_RBS_SS_0_8
from replace_agr import RBS_ONLY
from replace_agr import SS_ONLY
from replace_agr import TEST_RECODING__NO_CONSERVATION_SINGLE_RULE_DIR


PWD = os.path.dirname(os.path.realpath(__file__))

DATA_DIR = os.path.join(PWD, '../data/')

ORIGINAL_GENOME_DEST = os.path.join(PWD,
        '../data/reference_genomes/mds42_AP012306.1.genbank')

RECODED_GENOME_PATH = os.path.join(PWD,
        '../data/recoded_genomes/r_mds42_ap012306_agr_to_cgt.genbank')

MG1655_GENOME_PATH = os.path.join(PWD,
        '../data/reference_genomes/mg1655_NC_000913.3.genbank')

AGR_TO_ATG_DISTANCE = 21

AGR_CODONS = set(['AGG', 'AGA'])

NON_AGR_SYN_CODONS = set(['CGA', 'CGC', 'CGG', 'CGT'])

ARG_CODON_SET = AGR_CODONS | NON_AGR_SYN_CODONS

DOWNSTREAM_RBS_REPORT_FILE = os.path.join(PWD,
        '../data/analysis/genes_with_agr_in_rbs_of_downstream_rbs.csv')

AGR_IN_NON_ESSENTIAL_UPSTREAM_OF_ESSENTIAL = os.path.join(PWD,
        '../data/analysis/agr_in_non_essential_upstream_of_essential.csv')

RBS_IN_FIRST_30_NT_REPORT_FILE = os.path.join(PWD,
        '../data/analysis/genes_with_agr_in_first_30_nt.csv')

RBS_IN_FIRST_30_NT_REPORT_FILE__FORCE_FAKE_ATG = os.path.join(PWD,
        '../data/analysis/genes_with_agr_in_first_30_nt_force_fake_atg.csv')

CODON_COMPARISON_FIRST_30_NT = os.path.join(PWD,
        '../data/analysis/codon_comparison_AGR_first_30_nt.csv')

CODON_COMPARISON_RBS_OF_DOWNSTREAM = os.path.join(PWD,
        '../data/analysis/codon_comparison_AGR_in_RBS_of_downstream.csv')

CODON_COMPARISON_REMAINING = os.path.join(PWD,
        '../data/analysis/codon_comparison_remaining.csv')

CODON_COMPARISON_REMAINING_DEDUPED = os.path.join(PWD,
        '../data/analysis/codon_comparison_remaining_deduped.csv')

# NOTE: From getk_validation. Uses mix of Keio and PEC data.
KEIO_DATA = os.path.join(PWD,
    '../data/analysis/essential_genes.csv')


def _get_mg1655_gene_to_function_map():
    """The new mg1655 has awesome function annotations.

    We want to use these in our mds42 reports.
    """
    mg1655_seq_record = get_genome_record(MG1655_GENOME_PATH)
    coding_features = [f for f in mg1655_seq_record.features if f.type == 'CDS']

    gene_name_to_function_map = {}

    for f in coding_features:
        gene = get_feature_gene_name(f)
        function = get_feature_qualifier_value(f, 'function', '')
        gene_name_to_function_map[gene] = function

    return gene_name_to_function_map


def _get_last_AGR_base_position_in_feature(feature, seq_record):
    """Returns position relative to start of feature, 0-indexed.
    """
    f_seq = str(feature.extract(seq_record).seq).upper()
    assert len(f_seq) % 3 == 0
    last_AGR_base_position = -1
    for bp in range(0, len(f_seq), 3):
        codon = f_seq[bp:bp + 3]
        if codon in AGR_CODONS:
            last_AGR_base_position = bp
    return last_AGR_base_position


def _get_expression_comparison_data(feature, orig_seq, recoded_seq):
    rbs_calc_result_orig_seq = calc_rbs_score_for_feature(feature, orig_seq)
    rbs_calc_result_recoded_seq = calc_rbs_score_for_feature(
            feature, recoded_seq)

    expression_original = rbs_calc_result_orig_seq['expression']
    expression_recoded = rbs_calc_result_recoded_seq['expression']

    if expression_recoded == UNKNOWN_RBS_CALC_RESULT:
        delta_expression = INFINITY
        ratio = INFINITY
    else:
        # Calculate positive ratio.
        delta_expression = expression_recoded / expression_original
        if delta_expression >= 1:
            ratio = delta_expression
        else:
            ratio = 1.0 / delta_expression

    ratio = 1 + math.log(ratio)

    return {
        'expression_original': expression_original,
        'expression_recoded': expression_recoded,
        'delta_expression': delta_expression,
        'positive_ratio': ratio
    }


INFINITY = 99999999


def _insert_fake_ATG(seq, AGR_start, buffer_after_AGR=9):
    ATG_start = AGR_start + 3 + buffer_after_AGR
    new_seq = (
        seq[:ATG_start] +
        'ATG'
        # seq[ATG_start + 3:]
    )
    return new_seq


def _get_internal_rbs_expression_comparison_data(feature, orig_seq, recoded_seq,
        start, force_fake_ATG=False, buffer=20):
    orig_subseq = feature.extract(orig_seq)[2:start + buffer]
    if force_fake_ATG:
        orig_subseq = _insert_fake_ATG(orig_subseq, start - 2)
    rbs_calc_result_orig_seq = calc_internal_rbs_expression(orig_subseq,
            start - 2 + 3)

    recoded_subseq = feature.extract(recoded_seq)[2:start + buffer]
    if force_fake_ATG:
        recoded_subseq = _insert_fake_ATG(recoded_subseq, start - 2)
    rbs_calc_result_recoded_seq = calc_internal_rbs_expression(recoded_subseq,
            start - 2 + 3)

    expression_original = rbs_calc_result_orig_seq['expression']
    expression_recoded = rbs_calc_result_recoded_seq['expression']
    putative_start_orig = rbs_calc_result_orig_seq['putative_start']
    putative_start_recoded = rbs_calc_result_recoded_seq['putative_start']

    if (expression_original == UNKNOWN_RBS_CALC_RESULT and
            expression_recoded == UNKNOWN_RBS_CALC_RESULT):
        delta_expression = 1
    elif expression_original == UNKNOWN_RBS_CALC_RESULT:
        delta_expression = INFINITY
    elif expression_recoded == UNKNOWN_RBS_CALC_RESULT:
        delta_expression = 0
    elif expression_original == 0:
        if expression_recoded == 0:
            delta_expression = 1
        else:
            delta_expression = INFINITY
    else:
        delta_expression = expression_recoded / expression_original

    # Format subseqs for display.
    orig_subseq = _highlight_AGR(orig_subseq, start)
    recoded_subseq = _highlight_AGR(recoded_subseq, start)

    return {
        'rbs_calc_subseq_orig': orig_subseq,
        'rbs_calc_subseq_recoded': recoded_subseq,
        'expression_original': expression_original,
        'expression_recoded': expression_recoded,
        'delta_expression': delta_expression,
        'putative_start_orig': putative_start_orig,
        'putative_start_recoded': putative_start_recoded,
    }


def _highlight_AGR(orig_subseq, start):
    new_seq = (
        orig_subseq[:start - 2].lower() +
        orig_subseq[start - 2:start + 1].upper() +
        orig_subseq[start + 1:].lower()
    )
    assert new_seq.upper() == orig_subseq.upper()
    return new_seq


def analyze_genes_with_agr_in_rbs_of_downstream(seq_record, recoded_seq_record):
    orig_seq = str(seq_record.seq)
    recoded_seq = str(recoded_seq_record.seq)

    # Get mg1655 gene to function map so we can annotate function.
    gene_to_function_map = _get_mg1655_gene_to_function_map()

    # Use only coding features.
    coding_features = [f for f in seq_record.features if f.type == 'CDS']

    # Sort features by starting location.
    coding_features = sorted(coding_features, key=lambda f: f.location.start)

    # Keep track of flagged features.
    flagged_features = []

    # Identify flagged features. We only care about those facing in same
    # direction so we can handle them separately.
    fwd_coding_features = [f for f in coding_features if f.strand == 1]
    rev_coding_features = [f for f in coding_features if f.strand == -1]
    assert len(coding_features) == (
            len(fwd_coding_features) + len(rev_coding_features))

    # Check forward features.
    total_fwd_features = len(fwd_coding_features)
    for f_idx in range(total_fwd_features):
        if f_idx % 100 == 0:
            print "Analyzing forward feature %d of %d" % (f_idx, total_fwd_features)
        f = fwd_coding_features[f_idx]
        if f_idx >= len(fwd_coding_features) - 1:
            break
        next_f = fwd_coding_features[f_idx + 1]

        # Optimization to avoid AGR seek.
        if next_f.location.start - f.location.end > AGR_TO_ATG_DISTANCE:
            continue

        last_AGR_in_upstream_feature = _get_last_AGR_base_position_in_feature(
                f, seq_record)
        last_AGR_global_pos = last_AGR_in_upstream_feature + f.location.start
        if next_f.location.start - last_AGR_global_pos < AGR_TO_ATG_DISTANCE:
            flagged_locus_data = {
                'upstream_feature': f,
                'downstream_feature': next_f,
                'upstream_AGR_start_pos': last_AGR_in_upstream_feature + 1,
            }
            flagged_locus_data.update(_get_expression_comparison_data(
                    next_f, orig_seq, recoded_seq))
            flagged_features.append(flagged_locus_data)

    # Check reverse features.
    total_rev_features = len(rev_coding_features)
    for idx in range(total_rev_features - 1, 0, -1):
        if idx % 100 == 0:
            print "Analyzing rev feature %d of %d" % (idx, total_rev_features)
        f = rev_coding_features[idx]
        next_f = rev_coding_features[idx - 1]

        # Optimization to avoid AGR seek.
        if f.location.start - next_f.location.end > AGR_TO_ATG_DISTANCE:
            continue

        last_AGR_in_upstream_feature = _get_last_AGR_base_position_in_feature(
                f, seq_record)
        last_AGR_global_pos = (f.location.start + (
                len(f) - last_AGR_in_upstream_feature))
        if last_AGR_global_pos - next_f.location.end < AGR_TO_ATG_DISTANCE:
            flagged_locus_data = {
                'upstream_feature': f,
                'downstream_feature': next_f,
                'upstream_AGR_start_pos': last_AGR_in_upstream_feature + 1,
            }
            flagged_locus_data.update(_get_expression_comparison_data(
                    next_f, orig_seq, recoded_seq))
            flagged_features.append(flagged_locus_data)

    REPORT_FIELD_NAMES = [
        'upstream_gene_name',
        'upstream_gene_strand',
        'upstream_AGR_start_pos',
        'upstream_gene_function',
        'downstream_gene_name',
        'downstream_gene_strand',
        'downstream_gene_function',
        'expression_original',
        'expression_recoded',
        'delta_expression',
    ]

    with open(DOWNSTREAM_RBS_REPORT_FILE, 'w') as output_fh:
        writer = csv.DictWriter(output_fh, REPORT_FIELD_NAMES)
        writer.writeheader()
        for obj in flagged_features:
            uf = obj['upstream_feature']
            df = obj['downstream_feature']
            uf_gene = get_feature_gene_name(uf)
            df_gene = get_feature_gene_name(df)
            writer.writerow({
                'upstream_gene_name': uf_gene,
                'upstream_gene_strand': uf.strand,
                'upstream_AGR_start_pos': obj['upstream_AGR_start_pos'],
                'upstream_gene_function': gene_to_function_map.get(uf_gene, ''),
                'downstream_gene_name': df_gene,
                'downstream_gene_strand': df.strand,
                'downstream_gene_function': gene_to_function_map.get(df_gene, ''),
                'expression_original': obj['expression_original'],
                'expression_recoded': obj['expression_recoded'],
                'delta_expression': obj['delta_expression'],
            })


def _is_agr_in_first_30_nt(feature, seq_record):
    f_seq = str(feature.extract(seq_record).seq).upper()
    for bp in range(0, 30, 3):
        codon = f_seq[bp:bp + 3]
        if codon in AGR_CODONS:
            return bp
    return False


def _get_agr_pos_in_first_30_nt(feature, seq_record):
    agr_pos_list = []
    f_seq = str(feature.extract(seq_record).seq).upper()
    for bp in range(0, 30, 3):
        codon = f_seq[bp:bp + 3]
        if codon in AGR_CODONS:
            agr_pos_list.append((bp, codon))
    return agr_pos_list


def _get_all_agr_pos(feature, seq_record):
    agr_pos_list = []
    f_seq = str(feature.extract(seq_record).seq).upper()
    bp = 0
    while bp < len(feature):
        # HACK: Support for prfB frame shift.
        if get_feature_gene_name(feature) == 'prfB' and bp == 75:
            bp += 1
        codon = f_seq[bp:bp + 3]
        if codon in AGR_CODONS:
            agr_pos_list.append((bp, codon))
        bp += 3
    return agr_pos_list


def analyze_genes_with_AGR_in_first_30_nt(orig_seq_record, recoded_seq_record,
        output_report, force_fake_ATG=False):
    """Analyze genes with AGR in first 30 nt, as these might be significant
    in terms of pausing and/or unannotated RBS role.
    """
    orig_seq = str(orig_seq_record.seq).upper()
    recoded_seq = str(recoded_seq_record.seq).upper()

    # Get mg1655 gene to function map so we can annotate function.
    gene_to_function_map = _get_mg1655_gene_to_function_map()

    # Use only coding features.
    coding_features = [f for f in orig_seq_record.features if f.type == 'CDS']
    total_coding_features = len(coding_features)

    flagged_features = []

    for idx, f in enumerate(coding_features):
        if idx % 100 == 0:
            print "Analyzing feature %d of %d" % (idx, total_coding_features)
        agr_bp = _is_agr_in_first_30_nt(f, orig_seq_record)
        if agr_bp:
            flagged_locus_data = {
                'feature': f,
                'AGR_start_pos': agr_bp + 1
            }
            flagged_locus_data.update(
                    _get_internal_rbs_expression_comparison_data(
                            f, orig_seq, recoded_seq, agr_bp,
                            force_fake_ATG=True))
            flagged_features.append(flagged_locus_data)

    REPORT_FIELD_NAMES = [
        'gene',
        'strand',
        'AGR_start_pos',
        'function',
        'expression_original',
        'expression_recoded',
        'delta_expression',
        'rbs_calc_subseq_orig',
        'rbs_calc_subseq_recoded',
        'putative_start_orig',
        'putative_start_recoded'
    ]

    with open(output_report, 'w') as output_fh:
        writer = csv.DictWriter(output_fh, REPORT_FIELD_NAMES)
        writer.writeheader()
        for obj in flagged_features:
            f = obj['feature']
            gene = get_feature_gene_name(f)
            writer.writerow({
                'gene': gene,
                'strand': f.strand,
                'function': gene_to_function_map.get(gene, ''),
                'AGR_start_pos': obj['AGR_start_pos'],
                'expression_original': obj['expression_original'],
                'expression_recoded': obj['expression_recoded'],
                'delta_expression': obj['delta_expression'],
                'rbs_calc_subseq_orig': obj['rbs_calc_subseq_orig'],
                'putative_start_orig': obj['putative_start_orig'],
                'rbs_calc_subseq_recoded': obj['rbs_calc_subseq_recoded'],
                'putative_start_recoded': obj['putative_start_recoded']
            })


def generate_alternate_codon_choice_comparison_first_30_nt(
        seq_record, output_report, all_codons=False):
    """Generates scores for codon comparisons.

    Args:
        all_codons: If True, do more than just first 30nt.
    """
    data_list = []

    # Reusable scorers.
    refactor_context = _build_refactor_context(seq_record)
    ss_scorer = SSScorerOutsideFeature(refactor_context)
    rbs_scorer = PreserveRBSScorer(refactor_context)

    coding_features = [f for f in seq_record.features if f.type == 'CDS']
    total_coding_features = len(coding_features)

    # Analyze genes with AGR in first 30 bp.
    for idx, f in enumerate(coding_features):
        if idx % 100 == 0:
            print "Analyzing feature %d of %d" % (idx, total_coding_features)

        gene = get_feature_gene_name(f)

        if all_codons:
            find_AGR_fn = _get_all_agr_pos
        else:
            find_AGR_fn = _get_agr_pos_in_first_30_nt

        # AGR in first 30 bp.
        agr_bp_list = find_AGR_fn(f, seq_record)
        if agr_bp_list:
            for agr_pos, wt_codon in agr_bp_list:
                if agr_pos <= 30:
                    category = 'AGR in first 30'
                else:
                    category = 'remaining AGR'
                data_obj = {
                    'gene': gene,
                    'strand': f.location.strand,
                    'AGR_pos': agr_pos + 1,
                    'wt_codon': wt_codon,
                    'type': category
                }

                data_obj.update(
                        _get_profile_scores_at_agr_pos(
                                agr_pos, seq_record, f, ss_scorer, rbs_scorer))

                data_list.append(data_obj)

    # Write results.
    df = pd.DataFrame(data_list)
    column_order = [
        'gene',
        'strand',
        'AGR_pos',
        'wt_codon',
        'type'
    ]
    for codon in sorted(NON_AGR_SYN_CODONS):
        column_order.append('mRNA_' + codon)
    for codon in sorted(NON_AGR_SYN_CODONS):
        column_order.append('RBS_' + codon)
    df = df[column_order]
    df.to_csv(output_report, index=False)


def _get_profile_scores_at_agr_pos(
        agr_start_pos, seq_record, feature, ss_scorer, rbs_scorer):
    """Returns comparison dict for codon index.
    """
    data_dict = {}
    for codon in NON_AGR_SYN_CODONS:
        data_dict.update(
                _calculate_scores_for_alt_codon(
                        codon, agr_start_pos, seq_record, feature, ss_scorer,
                        rbs_scorer))
    return data_dict


def _calculate_scores_for_alt_codon(
        codon, agr_start_pos, seq_record, feature, ss_scorer, rbs_scorer):
    """Returns dictionary with profile scores for codon.
    """
    data_dict = {}

    BUFFER_CUT = 100

    orig_seq = str(seq_record.seq[
            feature.location.start - BUFFER_CUT:
            feature.location.end + BUFFER_CUT])

    # Determine codon start position relative to truncated seq.
    if feature.strand == 1:
        codon_start_idx = agr_start_pos + BUFFER_CUT
        assert orig_seq[codon_start_idx:codon_start_idx + 3] in ['AGG', 'AGA']

        if agr_start_pos <= 30:
            start_codon_idx = BUFFER_CUT
        else:
            start_codon_idx = codon_start_idx
    else:
        codon_start_idx = len(feature) - agr_start_pos + BUFFER_CUT

        # HACK: prfB
        if get_feature_gene_name(feature) == 'prfB' and agr_start_pos < 75:
            codon_start_idx += 1

        is_correct = (reverse_complement(orig_seq[
                codon_start_idx-3:codon_start_idx]) in ['AGG', 'AGA'])
        assert is_correct
        if agr_start_pos <= 30:
            start_codon_idx = len(orig_seq) - BUFFER_CUT
        else:
            start_codon_idx = codon_start_idx

    # Update the sequence.
    swap_codon = codon
    if feature.strand == 1:
        swap_codon = codon
    else:
        swap_codon = reverse_complement(codon)

    if feature.strand == 1:
        new_seq = (
                orig_seq[:codon_start_idx] +
                swap_codon +
                orig_seq[codon_start_idx + 3:]
        )
    else:
        new_seq = (
                orig_seq[:codon_start_idx - 3] +
                swap_codon +
                orig_seq[codon_start_idx:]
        )

    # Calculate mRNA score.
    ss_score = ss_scorer.score_at_codon_start_index(
            start_codon_idx, orig_seq, new_seq, feature, raw_ratio=True)
    ss_score = float("{0:.4f}".format(ss_score))
    data_dict['mRNA_' + codon] = ss_score

    # Calculate RBS score.
    if feature.strand == 1:
        polarity_adjust_orig_seq = orig_seq
        polarity_adjusted_new_seq = new_seq
        polarity_adjusted_codon_start_idx = codon_start_idx
    else:
        polarity_adjust_orig_seq = reverse_complement(orig_seq)
        polarity_adjusted_new_seq = reverse_complement(new_seq)
        polarity_adjusted_codon_start_idx = len(orig_seq) - codon_start_idx
    rbs_score_metadata = rbs_scorer.score_at_codon_start_index(
            polarity_adjusted_codon_start_idx, polarity_adjust_orig_seq,
            polarity_adjusted_new_seq, feature, result_with_metadata=True,
            buffer_after_AGR_interval_list=[12, 13])
    rbs_score = rbs_score_metadata['RBS_delta']
    rbs_score = float("{0:.4f}".format(rbs_score))
    data_dict['RBS_' + codon] = rbs_score

    return data_dict


def generate_alternate_codon_choice_comparison_rbs_of_downstream(
        seq_record, output_report):
    data_list = []

    # Use only coding features.
    coding_features = [f for f in seq_record.features if f.type == 'CDS']

    # Sort features by starting location.
    coding_features = sorted(coding_features, key=lambda f: f.location.start)

    fwd_coding_features = [f for f in coding_features if f.strand == 1]
    rev_coding_features = [f for f in coding_features if f.strand == -1]
    assert len(coding_features) == (
            len(fwd_coding_features) + len(rev_coding_features))

    # Copy of seq record for modifying.
    mod_seq_record = copy.deepcopy(seq_record)

    # Reusable scorers.
    refactor_context = _build_refactor_context(seq_record)
    ss_scorer = SSScorerOutsideFeature(refactor_context)

    # Check forward features.
    print 'Analyzing forward ...'
    total_fwd_features = len(fwd_coding_features)
    for f_idx in range(total_fwd_features):
        if f_idx % 100 == 0:
            print "Analyzing forward feature %d of %d" % (f_idx, total_fwd_features)
        if f_idx >= len(fwd_coding_features) - 1:
            break
        f = fwd_coding_features[f_idx]
        next_f = fwd_coding_features[f_idx + 1]

        # Optimization to avoid AGR seek.
        if next_f.location.start - f.location.end > AGR_TO_ATG_DISTANCE:
            continue

        last_AGR_in_upstream_feature = _get_last_AGR_base_position_in_feature(
                f, seq_record)
        last_AGR_global_pos = last_AGR_in_upstream_feature + f.location.start
        if next_f.location.start - last_AGR_global_pos > AGR_TO_ATG_DISTANCE:
            continue

        wt_codon = str(
                seq_record.seq[last_AGR_global_pos:last_AGR_global_pos + 3])
        assert wt_codon in AGR_CODONS

        data_obj = {
            'gene': get_feature_gene_name(f),
            'strand': f.location.strand,
            'downstream_gene': get_feature_gene_name(next_f),
            'AGR_pos': last_AGR_in_upstream_feature + 1,
            'wt_codon': wt_codon,
            'type': 'AGR in RBS of downstream'
        }

        try:
            for codon in NON_AGR_SYN_CODONS:
                mod_seq_record.seq = (
                        seq_record.seq[:last_AGR_global_pos] +
                        codon +
                        seq_record.seq[last_AGR_global_pos + 3:])

                # Compute RBS score.
                expression_compare_data = _get_expression_comparison_data(
                        next_f, seq_record.seq, mod_seq_record.seq)
                data_obj['RBS_' + codon] = float("{0:.4f}".format(
                        expression_compare_data['delta_expression']))

                # Compute mRNA free energy change.
                ss_score = ss_scorer.score_at_codon_start_index(
                        next_f.location.start, seq_record.seq,
                        mod_seq_record.seq, next_f, raw_ratio=True)
                data_obj['mRNA_' + codon] = float("{0:.4f}".format(ss_score))

        except TypeError:
            # Happens, for example, at frlC/frlD junction where changing AGR
            # in frlC breaks ATG of frlD.
            data_obj['RBS_' + codon] = INFINITY

        data_list.append(data_obj)

    # Check reverse features.
    print 'Analyzing reverse ...'
    total_rev_features = len(rev_coding_features)
    for idx in range(total_rev_features - 1, 0, -1):
        if idx % 100 == 0:
            print "Analyzing rev feature %d of %d" % (idx, total_rev_features)
        f = rev_coding_features[idx]
        next_f = rev_coding_features[idx - 1]

        # TODO: Special handling for multi-component genes.
        if get_feature_gene_name(f) == 'prfB':
            continue

        # Optimization to avoid AGR seek.
        if f.location.start - next_f.location.end > AGR_TO_ATG_DISTANCE:
            continue

        last_AGR_in_upstream_feature = _get_last_AGR_base_position_in_feature(
                f, seq_record)
        if last_AGR_in_upstream_feature == -1:
            continue
        last_AGR_global_pos = (f.location.start + (
                len(f) - last_AGR_in_upstream_feature))

        wt_codon = str(reverse_complement(
                seq_record.seq[last_AGR_global_pos - 3:last_AGR_global_pos]))
        assert wt_codon in AGR_CODONS
        data_obj = {
            'gene': get_feature_gene_name(f),
            'strand': f.location.strand,
            'downstream_gene': get_feature_gene_name(next_f),
            'AGR_pos': last_AGR_in_upstream_feature + 1,
            'wt_codon': wt_codon,
            'type': 'AGR in RBS of downstream'
        }

        for codon in NON_AGR_SYN_CODONS:
            mod_seq_record.seq = (
                    seq_record.seq[:last_AGR_global_pos - 3] +
                    reverse_complement(codon) +
                    seq_record.seq[last_AGR_global_pos:])

            # RBS score
            expression_compare_data = _get_expression_comparison_data(
                    next_f, seq_record.seq, mod_seq_record.seq)
            data_obj['RBS_' + codon] = float("{0:.4f}".format(
                    expression_compare_data['delta_expression']))

            # Compute mRNA free energy change.
            ss_score = ss_scorer.score_at_codon_start_index(
                    next_f.location.end, seq_record.seq,
                    mod_seq_record.seq, next_f, raw_ratio=True)
            data_obj['mRNA_' + codon] = float("{0:.4f}".format(ss_score))

        data_list.append(data_obj)

    # Write results.
    df = pd.DataFrame(data_list)
    column_order = [
        'gene',
        'strand',
        'downstream_gene',
        'AGR_pos',
        'wt_codon',
        'type'
    ]
    for codon in sorted(NON_AGR_SYN_CODONS):
        column_order.append('mRNA_' + codon)
    for codon in sorted(NON_AGR_SYN_CODONS):
        column_order.append('RBS_' + codon)
    df = df[column_order]
    df.to_csv(output_report, index=False)


def compare_to_designed_genomes(
        codon_comparison_record, orig_seq_record, recoded_genome_list,
        output_report, comparing_AGR_in_first_30_nt=True):
    """Follows up on codon comparison functions above.

    Args:
        codon_comparison_record: Path to file containing alternate codon scores.
        recoded_genome_list: List of paths to recoded genomes.
    """
    assert len(recoded_genome_list) == 1

    df = pd.read_csv(codon_comparison_record)


    # Copy of seq record for modifying.
    mod_seq_record = copy.deepcopy(orig_seq_record)

    best_mRNA_codons = []
    best_RBS_codons = []

    # Initialize dictionary with key genome name and value list of codon
    # choices.
    genome_dict = {}
    for genome_path in recoded_genome_list:
        genome_name = os.path.splitext(os.path.split(genome_path)[1])[0]
        seq_record = get_genome_record(genome_path)
        genome_dict[genome_name] = {
            'seq_record': seq_record,
            'codon_choice_list': [],
            'recoded_mRNA_score': [],
            'recoded_RBS_score': [],
        }

        # HACK: Previously we took a list but no longer. Anyway I am going
        # too quickly to clean this up or even explain it right now.
        # Reusable scorers.
        refactor_context = _build_refactor_context(seq_record)
        ss_scorer = SSScorerOutsideFeature(refactor_context)
        rbs_scorer = PreserveRBSScorer(refactor_context)

    assert len(genome_dict) == 1

    def _append_best_codon(row, prefix, growing_best_codon_list):
        best_mRNA_codon = 'CGG'
        best_mRNA_score = INFINITY
        all_equal = True
        last_val = None
        for codon in NON_AGR_SYN_CODONS:
            score = row[prefix + codon]

            # Rolling check to see if all equal.
            if last_val is not None:
                if score != last_val:
                    all_equal = False
            last_val = score

            # Check if best codon.
            if score < best_mRNA_score:
                best_mRNA_score = score
                best_mRNA_codon = codon

        # Set best codon.
        if all_equal:
            growing_best_codon_list.append('*')
        else:
            growing_best_codon_list.append(best_mRNA_codon)

    for idx, row in df.iterrows():

        # Select best mRNA.
        if 'mRNA_CGG' in df:
            _append_best_codon(row, 'mRNA_', best_mRNA_codons)

        # Select best RBS.
        _append_best_codon(row, 'RBS_', best_RBS_codons)

        # Select corresponding from each genome.
        for genome_name, genome_data in genome_dict.iteritems():
            # Identify which codon was chosen during recoding.
            f = get_feature_by_gene_name(
                    row['gene'], 'CDS', genome_data['seq_record'])
            seq_record = genome_data['seq_record']
            seq = genome_data['seq_record'].seq
            f_seq = f.extract(seq)
            AGR_pos = row['AGR_pos'] - 1
            codon = str(f_seq[AGR_pos:AGR_pos + 3])

            # Figure out scores for this codon inserted into original sequence.
            orig_feature = get_feature_by_gene_name(
                    row['gene'], 'CDS', orig_seq_record)

            # Different computation depending on gene context.
            if comparing_AGR_in_first_30_nt:
                score_dict = _calculate_scores_for_alt_codon(
                        codon, AGR_pos, orig_seq_record, orig_feature, ss_scorer,
                        rbs_scorer)
            else:
                score_dict = {}

                next_f = get_feature_by_gene_name(
                        row['downstream_gene'], 'CDS', orig_seq_record)

                # Modify sequence.

                if orig_feature.strand == 1:
                    last_AGR_global_pos = (
                            AGR_pos + orig_feature.location.start)
                    wt_codon = str(
                            orig_seq_record.seq[last_AGR_global_pos:
                                    last_AGR_global_pos + 3])
                    assert wt_codon in AGR_CODONS, wt_codon
                    swap_codon = codon
                else:
                    last_AGR_global_pos = (orig_feature.location.start + (
                        len(orig_feature) - AGR_pos))
                    rc_wt_codon = orig_seq_record.seq[
                            last_AGR_global_pos - 3:last_AGR_global_pos]
                    wt_codon = str(reverse_complement(rc_wt_codon))
                    assert wt_codon in AGR_CODONS, wt_codon
                    swap_codon = reverse_complement(codon)
                mod_seq_record.seq = (
                        orig_seq_record.seq[:last_AGR_global_pos] +
                        swap_codon +
                        orig_seq_record.seq[last_AGR_global_pos + 3:])

                # Compute RBS score.
                expression_compare_data = _get_expression_comparison_data(
                        next_f, orig_seq_record.seq, mod_seq_record.seq)
                score_dict['RBS_' + codon] = float("{0:.2f}".format(
                        expression_compare_data['positive_ratio']))

                # Compute mRNA free energy change.
                if next_f.strand == 1:
                    start_pos = next_f.location.start
                else:
                    start_pos = next_f.location.end
                ss_score = ss_scorer.score_at_codon_start_index(
                        start_pos, orig_seq_record.seq,
                        mod_seq_record.seq, next_f)
                score_dict['mRNA_' + codon] = float("{0:.2f}".format(ss_score))

            # DEBUG
            # score_dict = {
            #     'mRNA': 0,
            #     'RBS': 0,
            # }

            # DEBUG
            # print row['gene'], row['AGR_pos']
            # print (
            #     f_seq[AGR_pos - 9:AGR_pos].lower() +
            #     codon +
            #     f_seq[AGR_pos - 3:AGR_pos + 9].lower()
            # )
            # print codon

            # assert codon in NON_AGR_SYN_CODONS
            genome_data['codon_choice_list'].append(codon)

            for key, value in score_dict.iteritems():
                if re.match('mRNA', key):
                    genome_data['recoded_mRNA_score'].append(value)
                elif re.match('RBS', key):
                    genome_data['recoded_RBS_score'].append(value)
                else:
                    raise AssertionError("Unrecognized scores")

    if len(best_mRNA_codons):
        df['best_ss'] = best_mRNA_codons
    df['best_RBS'] = best_RBS_codons

    for genome_name, genome_data in genome_dict.iteritems():
        df[genome_name] = genome_data['codon_choice_list']
        df['recoded_mRNA_score'] = genome_data['recoded_mRNA_score']
        df['recoded_RBS_score'] = genome_data['recoded_RBS_score']

    df.to_csv(output_report, index=False)


def _analyze_agr_in_first_30_nt(g, orig_seq_record):
    ANALYSIS_AGR_FIRST_30NT = (
            os.path.splitext(g)[0] + '.analysis.AGR_first_30nt.csv')
    compare_to_designed_genomes(
            CODON_COMPARISON_FIRST_30_NT, orig_seq_record, [g],
            ANALYSIS_AGR_FIRST_30NT)


def _analyze_agr_in_rbs_of_downstream(g, orig_seq_record):
    ANALYSIS_AGR_IN_RBS_OF_DOWNSTREAM = (
            os.path.splitext(g)[0] +
            '.analysis.AGR_in_rbs_of_downstream.csv')
    compare_to_designed_genomes(
            CODON_COMPARISON_RBS_OF_DOWNSTREAM, orig_seq_record, [g],
            ANALYSIS_AGR_IN_RBS_OF_DOWNSTREAM,
            comparing_AGR_in_first_30_nt=False)


def _analyze_agr_replacement(g, original_seq_record):
    refactored_seq_record = get_genome_record(g)

    _check_recoding(original_seq_record, refactored_seq_record)

    output_report = (
            os.path.splitext(g)[0] +
            '.analysis.global_agr_report.csv')

    # Store a list of dictionary objects, one for each AGR. The dictionary
    # contains data for that AGR.
    data_list = []

    coding_features = [f for f in original_seq_record.features
        if f.type == 'CDS']
    total_coding_features = len(coding_features)

    # Analyze genes with AGR in first 30 bp.
    for idx, f in enumerate(coding_features):
        # Verbose output.
        if idx % 100 == 0:
            print "Analyzing feature %d of %d" % (idx, total_coding_features)

        new_feature = get_feature_by_id(refactored_seq_record, f.id)
        assert new_feature is not None
        new_seq = str(new_feature.extract(refactored_seq_record.seq))

        f_seq = str(f.extract(original_seq_record).seq).upper()
        for bp in range(0, len(f_seq), 3):
            codon = f_seq[bp:bp + 3]
            if codon in AGR_CODONS:

                data_obj = {
                    'gene': get_feature_gene_name(f),
                    'pos_in_gene': bp,
                    'orig_codon': codon,
                    'new_codon': new_seq[bp:bp + 3]
                    # 'strand': f.location.strand,
                    # 'AGR_pos': agr_pos + 1,
                    # 'type': 'AGR in first 30'
                }
                data_list.append(data_obj)

    df = pd.DataFrame(data_list)
    column_order = [
        'gene',
        'pos_in_gene',
        'orig_codon',
        'new_codon'
        # 'strand',
        # 'AGR_pos',
        # 'type'
    ]
    df = df[column_order]
    df.to_csv(output_report, index=False)


def remove_duplicates_from_codon_comparison_remaining():
    """The code to put together the codon comparison scores for codons that are
    neither in the first 30 nt nor in RBS of a downstream gene includes
    elements from other categories that may not be scored appropriately.
    Remove these.
    """
    other_df = pd.read_csv(CODON_COMPARISON_REMAINING)

    first_30_nt_df = pd.read_csv(CODON_COMPARISON_FIRST_30_NT)
    rbs_in_downstream_df = pd.read_csv(CODON_COMPARISON_RBS_OF_DOWNSTREAM)
    remove_dup_keys = (set(
        [tuple(x) for x in first_30_nt_df[['gene', 'AGR_pos']].values] +
        [tuple(x) for x in rbs_in_downstream_df[['gene', 'AGR_pos']].values]))

    def _is_not_dup(row):
        return not (row['gene'], row['AGR_pos']) in remove_dup_keys
    other_df = other_df[other_df.apply(_is_not_dup, axis=1)]

    other_df.to_csv(CODON_COMPARISON_REMAINING_DEDUPED, index=False)


def find_agr_in_non_essential_upstream_of_essential():
    agr_df = pd.read_csv(DOWNSTREAM_RBS_REPORT_FILE)
    keio_df = pd.read_csv(KEIO_DATA)

    essential_gene_set = set(keio_df['gene'])
    assert len(essential_gene_set) == 299

    def _is_agr_in_non_essential_upstream_of_essential(row):
        return (row['upstream_gene_name'] not in essential_gene_set and
                row['downstream_gene_name'] in essential_gene_set)

    agr_passing_filter_df = agr_df[agr_df.apply(
            _is_agr_in_non_essential_upstream_of_essential, axis=1)]

    agr_passing_filter_df.to_csv(
            AGR_IN_NON_ESSENTIAL_UPSTREAM_OF_ESSENTIAL, index=False)


def main():
    # original_seq_record = _get_original_genome_record()

    # generate_alternate_codon_choice_comparison_first_30_nt(
    #         original_seq_record, CODON_COMPARISON_FIRST_30_NT)

    # generate_alternate_codon_choice_comparison_rbs_of_downstream(
    #         original_seq_record, CODON_COMPARISON_RBS_OF_DOWNSTREAM)

    # generate_alternate_codon_choice_comparison_first_30_nt(
    #         original_seq_record, CODON_COMPARISON_REMAINING, all_codons=True)

    # remove_duplicates_from_codon_comparison_remaining()

    # _analyze_agr_in_first_30_nt(SS_ONLY)
    # _analyze_agr_in_rbs_of_downstream(SS_ONLY)

    # _analyze_agr_in_first_30_nt(RBS_ONLY)
    # _analyze_agr_in_rbs_of_downstream(RBS_ONLY)

    # ALL_MIX_RULES_GENOMES = [
    #         MIX_RBS_SS_0_8,
    #         MIX_RBS_SS_0_65,
    #         MIX_RBS_SS_0_2,
    #         MIX_RBS_SS_0_5
    # ]

    # for g in ALL_MIX_RULES_GENOMES:
    #     _analyze_agr_in_first_30_nt(g)
    #     _analyze_agr_in_rbs_of_downstream(g)

    # _analyze_agr_in_first_30_nt(
    #         MIX_NON_SYNONYMOUS_RBS_SS_0_65, original_seq_record)
    # _analyze_agr_in_rbs_of_downstream(
    #         MIX_NON_SYNONYMOUS_RBS_SS_0_65, original_seq_record)

    # Analyze final AGR usage.
    # _analyze_agr_replacement(
    #         MIX_RBS_SS_0_65, original_seq_record)

    find_agr_in_non_essential_upstream_of_essential()


if __name__ == '__main__':
    main()
