"""
Main script for replacing AGRs.
"""

from collections import defaultdict
import copy
import os
import re

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation

from getk.aggregate_scorer import ComboAggregateScorer
from getk.biopython_util import add_feature_to_seq_record
from getk.biopython_util import get_feature_by_id
from getk.biopython_util import get_feature_gene_name
from getk.biopython_util import get_feature_qualifier_value
from getk.biopython_util import get_genome_record
from getk.biopython_util import update_feature_seq_given_feature
from getk.codon_replacer import replace_codons_in_genome
from getk.codon_replacer import SimpleCodonReplacer
from getk.codon_usage_memex import get_ecoli_codon_usage_memex
from getk.iterators import CodonIterator
from getk.refactor_checker import check_all
from getk.refactor_checker import check_forbidden_codons_removed
from getk.refactor_checker import RNA_TYPES
from getk.refactor_context import RefactorContext
from getk.scorer import GCScorer
from getk.scorer import LengthScorer
from getk.scorer import HomologyScorer
from getk.scorer import MaintainTranslationScorer
from getk.scorer import PreserveRBSScorer
from getk.scorer import SSScorer
from getk.scorer import SSScorerOutsideFeature

from conservation import find_agr_intervals
from conservation import replace_AGR_in_genome_using_conservation


### File locations

PWD = os.path.dirname(os.path.realpath(__file__))

DATA_DIR = os.path.join(PWD, '../data/')

RECODED_GENOMES_DIR = os.path.join(DATA_DIR, 'recoded_genomes')

VERSIONS_SUBDIR = os.path.join(RECODED_GENOMES_DIR, 'agr_recoded_versions')

ORIGINAL_GENOME_PATH = os.path.join(DATA_DIR,
        'reference_genomes/mg1655_U00096.3_2014_08_01.gb')

RECODED_GENOME_PATH = os.path.join(RECODED_GENOMES_DIR,
        'agr_recoded.gb')

INT1_CONSERVATION_REPLACEMENTS_PATH = os.path.join(RECODED_GENOMES_DIR,
        'agr_recoded.int1_conservation_replacements.gb')

RECODING_STATS_REPORT = os.path.join(RECODED_GENOMES_DIR,
        'r_mg1655_U00096.3_2014_08_01.recoding_stats.txt')

RBS_PROFILE_CACHE_LOCATION = os.path.join(DATA_DIR, 'cache',
        'mg1655_U00096.3_2014_08_01.rbs_strength_profile.pickle')

FIX_OVERLAPS_CACHE_LOCATION = os.path.join(DATA_DIR, 'cache', 'fix_overlaps')

TESTING_RULES_DIR = os.path.join(RECODED_GENOMES_DIR, 'testing_rules')

TEST_RECODING__NO_CONSERVATION_SINGLE_RULE_DIR = os.path.join(
            TESTING_RULES_DIR, 'no_conservation_single_rule')

TEST_RECODING__NO_CONSERVATION_MIXED_RULES_DIR = os.path.join(
            TESTING_RULES_DIR, 'no_conservation_mixed_rules')

TEST_RECODING__NO_CONSERVATION__NON_SYNONYMOUS__MIXED_DIR = os.path.join(
            TESTING_RULES_DIR, 'no_conservation_non_synonymous_mixed_rules')

TEST_RECODING__NO_CONSERVATION__NON_SYNONYMOUS__MIXED_DIR__SEGMENT_CACHE = (
        os.path.join(TEST_RECODING__NO_CONSERVATION__NON_SYNONYMOUS__MIXED_DIR,
                'segment_cache'))

### Recoded genomes

RBS_ONLY = os.path.join(
        TEST_RECODING__NO_CONSERVATION_SINGLE_RULE_DIR,
        'test_rules.rbs.2015_06_11_1530.gb')

SS_ONLY = os.path.join(
        TEST_RECODING__NO_CONSERVATION_SINGLE_RULE_DIR,
        'test_rules.ss_minus40_plus90.2015_06_10_1415.gb')

MIX_RBS_SS_0_8 = os.path.join(
        TEST_RECODING__NO_CONSERVATION_MIXED_RULES_DIR,
        'test_rules.ss_0.8_rbs_0.2.2015_06_10_1508.gb')

MIX_RBS_SS_0_65 = os.path.join(
        TEST_RECODING__NO_CONSERVATION_MIXED_RULES_DIR,
        'test_rules.ss_0.65_over_1.5411_rbs_0.35_over_8.4257.2015_06_10_1545.gb')

MIX_RBS_SS_0_2 = os.path.join(
        TEST_RECODING__NO_CONSERVATION_MIXED_RULES_DIR,
        'test_rules.ss_0.2_rbs_0.8.2015_06_10_1546.gb')

MIX_RBS_SS_0_5 = os.path.join(
        TEST_RECODING__NO_CONSERVATION_MIXED_RULES_DIR,
        'test_rules.ss_0.5_rbs_0.5.2015_06_10_1547.gb')

MIX_NON_SYNONYMOUS_RBS_SS_0_65 = os.path.join(
        TEST_RECODING__NO_CONSERVATION__NON_SYNONYMOUS__MIXED_DIR,
        'test_rules.non_synonymous.ss_0.65_over_1.5411_rbs_0.35_over_8.4257.2015_06_14_2230.gb')


def _generate_version_genbank_filename(version):
    return os.path.join(VERSIONS_SUBDIR,
            'agr_recoded_v{version}.gb'.format(version=version))


### Functions

# TODO: Actually check. Right now just skipping.
# Genes that need special handling during refactor checking.
REFACTOR_CHECKER_SPECIAL_CASES = {
    # Internal stop codon.
    'fdhF': {
        'type': 'temp_swap',
        # 'position': 585,
        # 'from': 'TGA',
        # 'to': 'GGG'
    },

    # Internal stop codon.
    'fdnG': {
        'type': 'temp_swap',
        # 'position': 585,
        # 'from': 'TGA',
        # 'to': 'GGG'
    },

    # Internal stop codon.
    'fdoG': {
        'type': 'temp_swap',
        # 'position': 585,
        # 'from': 'TGA',
        # 'to': 'GGG'
    }
}

# Codons we are removing.
AGR_CODONS = set(['AGA', 'AGG'])

# List of gene prefixes to skip.
SKIP_GENES_PREFIX_LIST = [
    'ins'
]


def skip_feature_if_matches_skip_prefix(feature):
    """Returns True if feature gene name matches prefix.
    """
    feature_gene_name = get_feature_gene_name(feature)
    if feature_gene_name is None:
        return False
    for g_prefix in SKIP_GENES_PREFIX_LIST:
        if re.match(g_prefix, feature_gene_name):
            return True
    return False


def skip_feature_if_pseudogene(feature):
    """Returns True if feature gene name matches prefix.
    """
    feature_note = get_feature_qualifier_value(feature, 'note')
    if feature_note is None:
        return False
    return 'pseudogene' in feature_note


def skip_feature_if_prophage(feature):
    """Returns True if feature gene name matches prefix.
    """
    feature_note = get_feature_qualifier_value(feature, 'product')
    if feature_note is not None and 'prophage' in feature_note:
        return True

    feature_note = get_feature_qualifier_value(feature, 'note')
    if feature_note is None:
        return False
    return 'prophage' in feature_note


def _get_feature_seq_with_AGR_replaced_by_CGT(seq_feature, seq_record):
    """Returns string sequence for the feature with AGR's replaced by CGT.

    Sequence is in 5'-to-3' direction, e.g. ATG...TAG
    """
    orig_seq = str(seq_feature.extract(seq_record.seq))
    new_seq = ''
    for codon_idx in range(0, len(orig_seq), 3):
        codon = orig_seq[codon_idx:codon_idx + 3]
        if codon in AGR_CODONS:
            new_codon = 'CGT'
        else:
            new_codon = codon
        new_seq += new_codon
    return new_seq


def _get_coding_features(seq_record):
    """Returns the coding features that we refactor.
    """
    coding_features = []
    for f in seq_record.features:
        # Run feature through gauntlet. Only append if make it all the way.
        if f.type != 'CDS':
            continue

        gene_name = get_feature_gene_name(f)
        assert gene_name, f
        if re.match('ins', gene_name):
            continue

        coding_features.append(f)
    return coding_features


def _naive_agr_to_cgt(original_seq_record, refactored_seq_record):
    """Replace all AGR to CGT.

    NOTES:
        * Fails on overlaps.
        * Currently unused.
    """
    coding_features = _get_coding_features(original_seq_record)

    # Temporary: Iterate through coding features and replace all AGR
    # occurrences with CTG.
    for f in coding_features:
        if len(f) % 3 != 0:
            print "Feature len for %s not multiple of 3." % (
                    get_feature_gene_name(original_seq_record, f))
            continue
        new_feature_seq = _get_feature_seq_with_AGR_replaced_by_CGT(
                f, original_seq_record)
        update_feature_seq_given_feature(refactored_seq_record, f,
                new_feature_seq)


def _refactor_with_getk(refactored_seq_record):
    """Uses latest GETK libraries to do recoding.

    Mutates refactored_seq_record.
    """
    refactor_context = _build_refactor_context(refactored_seq_record)

    # Construct iterators.
    iterator = CodonIterator(refactor_context)

    # Construct scorers.
    hard_scorer_list = [
            # MaintainTranslationScorer(refactor_context)
    ]
    # hard_scorer_coef_list = [1]
    hard_scorer_coef_list = []

    soft_scorer_list = [
            SSScorerOutsideFeature(refactor_context),
            PreserveRBSScorer(refactor_context),
    ]
    # soft_scorer_coef_list = [1]
    # soft_scorer_coef_list = [0.8, 0.2]
    soft_scorer_coef_list = [0.65/1.5411, 0.35/8.4257]
    # soft_scorer_coef_list = [0.2, 0.8]
    # soft_scorer_coef_list = [0.5, 0.5]

    aggregate_scorer = ComboAggregateScorer(
            hard_scorer_list, hard_scorer_coef_list,
            soft_scorer_list, soft_scorer_coef_list)

    refactor_context.set_aggregate_scorer(aggregate_scorer)

    return replace_codons_in_genome(
            refactor_context, iterator, aggregate_scorer,
            # start_index=0,
            # end_index=100000
    )


def _build_refactor_context(orig_seq_record):
    """Returns RefactorContext object.
    """
    codon_usage_memex = get_ecoli_codon_usage_memex(randomized=True)
    params = {
        'original_seq': orig_seq_record,
        'codon_usage_memex': codon_usage_memex,
        'forbidden_codons': AGR_CODONS,
        'rbs_strength_profile_cache_location': RBS_PROFILE_CACHE_LOCATION,
        # 'parallelized': True,
        # 'num_cores': 4,
        # 'cache_overlap_fixes': False,
        # 'fix_overlaps_cache_location': FIX_OVERLAPS_CACHE_LOCATION
        'allow_non_synonymous_swaps': True,
        'cache_recoded_segments': True,
        'cache_recoded_segments_location': TEST_RECODING__NO_CONSERVATION__NON_SYNONYMOUS__MIXED_DIR__SEGMENT_CACHE
    }
    return RefactorContext(params)


def _recoding_stats(original_seq_record, refactored_seq_record):
    """Simple function to compare genomes.
    """
    coding_features = [f for f in original_seq_record.features
            if f.type == 'CDS']
    codon_swap_map = defaultdict(lambda: defaultdict(lambda: 0))
    for f in coding_features:
        orig_seq = str(f.extract(original_seq_record.seq))
        new_feature = get_feature_by_id(refactored_seq_record, f.id)
        new_seq = str(new_feature.extract(refactored_seq_record.seq))
        for codon_idx in range(0, len(orig_seq), 3):
            codon = orig_seq[codon_idx:codon_idx + 3]
            new_codon = new_seq[codon_idx:codon_idx + 3]
            if codon != new_codon:
                codon_swap_map[codon][new_codon] += 1
    # for codon, new_codon_map in codon_swap_map.iteritems():
    #     print codon, new_codon_map
    return codon_swap_map



def _annotate_codon_changes(original_seq_record, refactored_seq_record):
    """Annotates non-AGR codon changes for analysis.
    """
    assert False
    AGR_CHANGE = 'agr_change'
    NON_AGR_CODON_CHANGE_TYPE = 'non_agr_change'

    refactor_context = _build_refactor_context(original_seq_record)

    # scorer_list = [
    #     # GCScorer(refactor_context),
    #     # PreserveRBSScorer(refactor_context)
    #     # SSScorer(refactor_context),
    #     # LengthScorer(refactor_context),
    #     # HomologyScorer(refactor_context)
    # ]

    coding_features = [f for f in original_seq_record.features
            if f.type == 'CDS']
    total_coding_features = len(coding_features)

    for idx, f in enumerate(coding_features):
        if idx % 100 == 0:
            print "Analyzing feature %d of %d" % (idx, total_coding_features)
        orig_seq = str(f.extract(original_seq_record.seq))
        new_feature = get_feature_by_id(refactored_seq_record, f.id)
        new_seq = str(new_feature.extract(refactored_seq_record.seq))
        for codon_idx in range(0, len(orig_seq), 3):
            codon = orig_seq[codon_idx:codon_idx + 3]
            new_codon = new_seq[codon_idx:codon_idx + 3]
            if codon != new_codon:
                if f.strand == 1:
                    start = new_feature.location.start + codon_idx
                else:
                    start = new_feature.location.start + (
                            len(f) - codon_idx - 3)
                if codon in AGR_CODONS:
                    new_f_type = AGR_CHANGE
                else:
                    new_f_type = NON_AGR_CODON_CHANGE_TYPE
                new_f = SeqFeature(
                        location=FeatureLocation(start, start + 3),
                        type=new_f_type,
                        strand=f.strand,
                        id=NON_AGR_CODON_CHANGE_TYPE + '_' + str(start)
                )

                new_f.qualifiers['orig_codon'] = codon
                new_f.qualifiers['new_codon'] = new_codon

                _annotate_alternate_codon_scores(
                        new_f, new_codon, codon_idx, orig_seq, new_seq,
                        refactor_context, scorer_list)

                # DEBUG
                # print get_feature_gene_name(f), codon_idx, codon, new_codon, non_agr_f.location.start, non_agr_f.location.end

                try:
                    add_feature_to_seq_record(refactored_seq_record, new_f)
                except AssertionError:
                    continue


def _annotate_alternate_codon_scores(new_f, new_codon, codon_start_idx, orig_seq,
        new_seq, refactor_context, scorer_list):
    """Adds annotations for alternate codon scores at given codon position.
    """
    syn_codons = refactor_context.codon_usage_memex.get_synonymous_codons(
            new_codon)
    for alt_codon in syn_codons:
        assert new_codon == new_seq[codon_start_idx:codon_start_idx + 3]
        mod_seq = (new_seq[:codon_start_idx] + alt_codon +
                new_seq[codon_start_idx + 3:])
        for scorer in scorer_list:
            key = 'score_{scorer_name}_{codon}'.format(
                    scorer_name=type(scorer).get_name(),
                    codon=alt_codon)
            # Modify the sequence to swap in alternate codon.
            score = scorer.score_at_codon_start_index(
                    codon_start_idx, orig_seq, mod_seq)
            if score == 0.0:
                score = "0.0"
            assert score is not None
            new_f.qualifiers[key] = score


def _get_original_genome_record():
    return get_genome_record(ORIGINAL_GENOME_PATH,
            ignore_features_fn_list=[
                    skip_feature_if_matches_skip_prefix,
                    skip_feature_if_pseudogene,
                    skip_feature_if_prophage],
    )

# Features that we know are broken during recoding but we have verified
# manually.
KNOWN_BROKEN_FEATURE_ID_LIST = [
    'ncRNA_b4425'
]


def _check_recoding(original_seq_record, refactored_seq_record):
    check_all(original_seq_record, refactored_seq_record,
            special_cases=REFACTOR_CHECKER_SPECIAL_CASES,
            ignore_problems_in_feature_ids=KNOWN_BROKEN_FEATURE_ID_LIST
    )
    check_forbidden_codons_removed(refactored_seq_record, AGR_CODONS)


def _count_AGR_occurrences(seq_record):
    cds_features = [f for f in seq_record.features if f.type == 'CDS']
    count = 0
    for f in cds_features:
        count += len(find_agr_intervals(f, seq_record))
    return count


def recode(output_path):
    original_seq_record = _get_original_genome_record()

    # Make a copy so that we can keep the original for comparison.
    # Downstream functions should not be mutating the original but we haven't
    # confirmed that.
    # refactored_seq_record = copy.deepcopy(original_seq_record)
    # print 'AGR count before', _count_AGR_occurrences(refactored_seq_record)
    refactored_seq_record = original_seq_record

    # # Use conservation data to swap codons manually.
    # replace_AGR_result = replace_AGR_in_genome_using_conservation(
    #         refactored_seq_record,
    #         ignore_translation_gene_list=REFACTOR_CHECKER_SPECIAL_CASES,
    #         verbose=True)
    # print 'total', replace_AGR_result['total_AGRs']
    # print 'found', replace_AGR_result['found']
    # refactored_seq_record = replace_AGR_result['seq_record']
    # print 'AGR count after', _count_AGR_occurrences(refactored_seq_record)
    # with open(INT1_CONSERVATION_REPLACEMENTS_PATH, 'w') as output_fh:
    #     SeqIO.write(refactored_seq_record, output_fh, 'genbank')
    # return

    # refactored_seq_record = get_genome_record(
    #         INT1_CONSERVATION_REPLACEMENTS_PATH)

    # # Make sure that change is actually made.
    # assert str(refactored_seq_record.seq) != str(original_seq_record.seq)

    # # Check constraints maintained.
    # check_all(original_seq_record, refactored_seq_record,
    #         special_cases=REFACTOR_CHECKER_SPECIAL_CASES,
    #         ignore_problems_in_feature_ids=KNOWN_BROKEN_FEATURE_ID_LIST
    # )

    # Continue with GETK recoding.
    refactor_result = _refactor_with_getk(refactored_seq_record)
    refactored_seq_record = refactor_result['new_seq']

    # Write the result genome.
    with open(output_path, 'w') as output_fh:
        SeqIO.write(refactored_seq_record, output_fh, 'genbank')

    # with open(RECODING_STATS_REPORT, 'w') as recoding_stats_fh:
    #     codon_swap_map = _recoding_stats(
    #             original_seq_record, refactored_seq_record)
    #     for codon, new_codon_map in codon_swap_map.iteritems():
    #         recoding_stats_fh.write(codon)
    #         recoding_stats_fh.write('\n')
    #         recoding_stats_fh.write(str(new_codon_map))
    #         recoding_stats_fh.write('\n>>>>>>>>>>>>>>\n')

    # Check recoding.
    _check_recoding(original_seq_record, refactored_seq_record)


def main():
    print 'Start ...'

    # Single-rule

    # recode(RBS_ONLY)
    # recode(SS_ONLY)

    # Mixed

    # recode(MIX_RBS_SS_0_8)
    # recode(MIX_RBS_SS_0_65)
    # recode(MIX_RBS_SS_0_2)
    # recode(MIX_RBS_SS_0_5)

    # Non-synonymous codons.
    recode(MIX_NON_SYNONYMOUS_RBS_SS_0_65)


def analyze():
    original_seq_record = _get_original_genome_record()
    refactored_seq_record = get_genome_record(RECODED_GENOME_PATH)

    _check_recoding(original_seq_record, refactored_seq_record)

    with open(RECODING_STATS_REPORT, 'w') as recoding_stats_fh:
        codon_swap_map = _recoding_stats(
                original_seq_record, refactored_seq_record)
        for codon, new_codon_map in codon_swap_map.iteritems():
            recoding_stats_fh.write(codon)
            recoding_stats_fh.write('\n')
            recoding_stats_fh.write(str(new_codon_map))
            recoding_stats_fh.write('\n>>>>>>>>>>>>>>\n')

    # # Compare RNAs
    # original_rna_features = [f.id for f in original_seq_record.features
    #         if f.type in RNA_TYPES]
    # refactored_rna_features = [f.id for f in refactored_seq_record.features
    #         if f.type in RNA_TYPES]
    # print set(original_rna_features) - set(refactored_rna_features)


def add_annotations_to_recoded_genome():
    """Run after recoding.
    """
    VERSION = 2.4

    original_seq_record = _get_original_genome_record()
    refactored_seq_record = get_genome_record(RECODED_GENOME_PATH)
    refactored_seq_record.name = 'agr_recoded_v{version}'.format(
            version=VERSION)

    # Profile
    import cProfile
    from datetime import datetime
    TMP_FILE_PREFIX = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    CPROFILE_OUTPUT_DEST = os.path.join('profiling',
            TMP_FILE_PREFIX + '_cprofile.out')
    cProfile.runctx(
            '_annotate_codon_changes(original_seq_record, refactored_seq_record)',
            globals(), locals(),
            CPROFILE_OUTPUT_DEST)
    # _annotate_codon_changes(original_seq_record, refactored_seq_record)

    # Clean up RBS-Calc.
    from getk.rbs_calc import NuPACK
    if os.path.exists(NuPACK.current_dir):
        import shutil
        shutil.rmtree(NuPACK.current_dir)

    with open(_generate_version_genbank_filename(VERSION), 'w') as output_fh:
        SeqIO.write(refactored_seq_record, output_fh, 'genbank')


def analyze_test_ss_conservative_split():
    """Analyzing running recoding with ss scoring only.
    """
    original_seq_record = _get_original_genome_record()

    SS_ONLY__LONG = os.path.join(RECODED_GENOMES_DIR, 'test_rules.ss_100bp.gb')
    refactored_seq_record = get_genome_record(SS_ONLY__LONG)

    original_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            original_seq_record.features)
    refactored_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            refactored_seq_record.features)

    original_cds_ids = [f.id for f in original_coding_features]
    refactored_cds_ids = [f.id for f in refactored_coding_features]

    print set(original_cds_ids) - set(refactored_cds_ids)


def fix_test_ss_conservative():
    ALL_BUT_LAST_PATH = os.path.join(
            RECODED_GENOMES_DIR, 'test_rules.ss_100bp.all_but_last.gb')
    all_but_last = get_genome_record(ALL_BUT_LAST_PATH)

    LAST_DIV_PATH = os.path.join(
            RECODED_GENOMES_DIR, 'test_rules.ss_100bp.last_div.gb')
    last_div = get_genome_record(LAST_DIV_PATH)

    combined = all_but_last + last_div

    # Write the result genome.
    SS_ONLY__LONG = os.path.join(RECODED_GENOMES_DIR, 'test_rules.ss_100bp.gb')
    with open(SS_ONLY__LONG, 'w') as output_fh:
        SeqIO.write(combined, output_fh, 'genbank')

    # Check recoding.
    original_seq_record = _get_original_genome_record()
    _check_recoding(original_seq_record, combined)


def generate_AGR_pos_list():
    original_seq_record = _get_original_genome_record()

    coding_features = [f for f in original_seq_record.features
            if f.type == 'CDS']

    with open('agr_pos_list.csv', 'w') as fh:
        for f in coding_features:
            orig_seq = str(f.extract(original_seq_record.seq))
            for codon_idx in range(0, len(orig_seq), 3):
                codon = orig_seq[codon_idx:codon_idx + 3]
                if codon in AGR_CODONS:
                    if f.strand == 1:
                        start = f.location.start + codon_idx
                    else:
                        start = f.location.start + (
                                len(f) - codon_idx - 3)
                    fh.write(str(start) + '\n')


if __name__ == '__main__':
    main()
    # analyze()
    # add_annotations_to_recoded_genome()

    # refactored_seq_record = get_genome_record(
    #         _generate_version_genbank_filename(2.3))
    # with open(_generate_version_genbank_filename(2.4), 'w') as output_fh:
    #     SeqIO.write(refactored_seq_record, output_fh, 'genbank')
    # analyze_test_ss_conservative_split()
    # fix_test_ss_conservative()
    # generate_AGR_pos_list()
