"""
Functions that allow us to use conservation data for recoding.
"""

import copy
import os
import re
from uuid import uuid4

from Bio.Seq import reverse_complement
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
import pyliftover

from getk.biopython_util import add_feature_to_seq_record
from getk.biopython_util import get_genome_record
from getk.biopython_util import update_feature_seq_given_feature
from getk.divider import Divider
from getk.refactor_checker import check_all

# Used for debug.
# class FakeException(Exception):
#     pass


### File locations

PWD = os.path.dirname(os.path.realpath(__file__))

CONSERVATION_DATA_DIR = os.path.join(PWD, '../data/conservation')

ECOLI_CONSERVATION_MAF_ORIGINAL = os.path.join(
        CONSERVATION_DATA_DIR, 'ecoli_conservation.maf')

ECOLI_CONSERVATION_MAF_INT_1 = os.path.join(
        CONSERVATION_DATA_DIR, 'ecoli_conservation.intermediate_1.maf')

ECOLI_CONSERVATION_MAF_UPDATED = os.path.join(
        CONSERVATION_DATA_DIR, 'ecoli_conservation.updated_u00096.3.maf')

MG1655_GENBANK = os.path.join(
        PWD, '../data/reference_genomes/mg1655_U00096.3_2014_08_01.gb')

MG1655_U00096_2_GENBANK = os.path.join(
        PWD, '../data/reference_genomes/mg1655_U00096.2.gb')

LIFTOVER_CHAIN_FILE = os.path.join(
        PWD, '../data/reference_genomes/liftover/',
        'mg1655_U00096.2.to.mg1655_U00096.3_2014_08_01.liftOver')


### Constants

MG1655_MAF_STRAIN_NAME = 'eschColi_K12.chr'

CONSERVATION_CODON_SWAP_TYPE = 'conservation'


### Functions

def update_mg1655_positions(input_maf_path, output_maf_path):
    """Update MG1655 positions in the conservation data.

    The positions from the microbes.ucsc.edu conservation track are relative
    to an older version of MG1655. We need to update these to be relative to
    the latest genome.

    Outputs an update version of the maf to the output path.
    """
    _update_using_liftover(input_maf_path, output_maf_path)


def _update_using_liftover(input_maf_path, output_maf_path):
    """Update using liftover.
    """
    liftover_obj = pyliftover.LiftOver(LIFTOVER_CHAIN_FILE)
    # First, remove duplicates.
    _remove_duplicate_alignments(input_maf_path, ECOLI_CONSERVATION_MAF_INT_1)

    # chunk_generator = _generate_maf_chunks(ECOLI_CONSERVATION_MAF_INT_1)
    chunk_generator = _generate_maf_chunks(ECOLI_CONSERVATION_MAF_INT_1)
    with open(output_maf_path, 'w') as output_fh:
        while True:
            try:
                chunk_data = chunk_generator.next()

                if chunk_data['type'] == 'HEADER':
                    output_fh.write(chunk_data['chunk'])
                    continue

                # Grab the k12 start.
                lines = chunk_data['chunk'].split('\n')
                parsed_k12_line = _parse_alignment_line(lines[1])
                assert parsed_k12_line['strain'] == 'eschColi_K12.chr'

                # Liftover the start position.
                # convert_coordinate() returns a list of tuples where the 2nd
                # el of each tuple is the converted position. We'll ignore
                # everything else for simplicity for now.
                start = _get_ecoli_position_from_maf_line(lines[1])
                liftover_pos_list = liftover_obj.convert_coordinate(
                        'U00096.2', start)
                if not liftover_pos_list:
                    print 'Mapping not found for %d' % start
                    continue
                updated_start = str(liftover_pos_list[0][1])

                # Update the chunk to write.
                parsed_k12_line['start'] = updated_start
                updated_k12_line = _join_alignment_line_parts(parsed_k12_line)
                lines[1] = updated_k12_line
                updated_chunk = '\n'.join(lines)

                # Write chunk.
                output_fh.write(updated_chunk)
            except StopIteration:
                break


def _parse_alignment_line(line):
    """Returns dictionary representing parsed line.
    """
    s, strain, start, size, sense, genome_size, seq = line.split()
    assert s == 's'
    return {
        'strain': strain,
        'start': int(start),
        'size': int(size),
        'sense': sense,
        'genome_size': int(genome_size),
        'seq': seq
    }


def _join_alignment_line_parts(parts):
    """Joins parts into whole line.
    """
    return '\t'.join([
            's',
            parts['strain'],
            str(parts['start']),
            str(parts['size']),
            str(parts['sense']),
            str(parts['genome_size']),
            parts['seq']])


def _get_ecoli_position_from_maf_line(line):
    """Helper that returns int position of ecoli line.
    """
    parsed_line = _parse_alignment_line(line)
    assert parsed_line['strain'] == 'eschColi_K12.chr', parsed_line['strain']
    return parsed_line['start']


def _remove_duplicate_alignments(input_maf_path, output_maf_path):
    """Removes duplicate alignment positions, picking the the best one.
    """

    def _get_score_from_line(line):
        """Helper that returns score.

        NOTE: Parses as int to make comparison safer.

        Raises AssertionError if not alignment line.
        """
        assert line[0] == 'a', line
        return int(float(re.search(r'[0-9\.]+', line).group()))

    # First go through and build a map from start position to highest
    # alignment score.
    map_start_to_alignment_score = {}
    with open(input_maf_path) as input_fh:

        # Use state machine strategy.
        class STATES:
            SEEK_NEXT_ALIGNMENT_SCORE = 1
            SEEK_K12_POSITION_FOR_SCORE = 2

        state = STATES.SEEK_NEXT_ALIGNMENT_SCORE
        current_alignment_score = None

        for line in input_fh:
            if state == STATES.SEEK_NEXT_ALIGNMENT_SCORE:
                if line[0] != 'a':
                    continue
                current_alignment_score = _get_score_from_line(line)
                state = STATES.SEEK_K12_POSITION_FOR_SCORE
            elif state == STATES.SEEK_K12_POSITION_FOR_SCORE:
                start_pos = _get_ecoli_position_from_maf_line(line)
                score = map_start_to_alignment_score.get(start_pos, -1)
                if current_alignment_score > score:
                    map_start_to_alignment_score[start_pos] = (
                            current_alignment_score)
                current_alignment_score = None
                state = STATES.SEEK_NEXT_ALIGNMENT_SCORE

    # Now go through and only keep the duplicate alignment with best score.
    def _maybe_write_previous_chunk(chunk, fh, written_position_set):
        """Helper function that determines whether to write chunk and then
        writes it if necessary.
        """
        # Check whether position/score is in the map, meaning it's the highest
        # alignment score and the one that we want.
        lines = chunk.split('\n')
        score = _get_score_from_line(lines[0])
        start_pos = _get_ecoli_position_from_maf_line(lines[1])
        # Write only if highest score and has not been written before.
        if (map_start_to_alignment_score[start_pos] == score and
                not start_pos in written_position_set):
            fh.write(current_alignment_chunk)
            written_position_set.add(start_pos)

    # Even though we are picking the highest score, occasionally we'll get
    # the same score more than once. As far as I can tell, these correspond
    # to the same entry. Use this structure to keep track of positions that
    # have already been written.
    written_position_set = set()

    with open(output_maf_path, 'w') as output_fh:

        # Use state machine strategy.
        class WRITING_STATES:
            WRITING_HEADER = 0
            WRITING_CONTENT = 1

        state = WRITING_STATES.WRITING_HEADER
        current_alignment_chunk = ''

        with open(input_maf_path) as input_fh:
            for line in input_fh:
                if state == WRITING_STATES.WRITING_HEADER and line[0] != 'a':
                    output_fh.write(line)
                elif line[0] == 'a':
                    # Permanently change state. Remaining lines come in chunks
                    # representing an alignment.
                    state = WRITING_STATES.WRITING_CONTENT

                    # Arrived at start of alignment chunk.
                    if current_alignment_chunk:
                        _maybe_write_previous_chunk(
                                current_alignment_chunk, output_fh,
                                written_position_set)
                    # Start new alignment chunk.
                    current_alignment_chunk = line
                else:
                    # Append line to accumulating chunk.
                    current_alignment_chunk += line

        # Handle last chunk.
        _maybe_write_previous_chunk(
                current_alignment_chunk, output_fh, written_position_set)

    print 'BEGIN ASSERT'

    # Assert that all positions written at least once.
    recorded_positions = set()

    chunk_generator = _generate_maf_chunks(output_maf_path)
    while True:
        try:
            chunk_data = chunk_generator.next()

            if chunk_data['type'] == 'HEADER':
                continue

            # Grab the k12 start.
            lines = chunk_data['chunk'].split('\n')
            parsed_k12_line = _parse_alignment_line(lines[1])
            assert parsed_k12_line['strain'] == 'eschColi_K12.chr'
            recorded_positions.add(parsed_k12_line['start'])
        except StopIteration:
            break

    expected_position_set = set(map_start_to_alignment_score.keys())
    print len(expected_position_set), len(recorded_positions)
    assert recorded_positions == expected_position_set, expected_position_set - recorded_positions


def _generate_maf_chunks(input_maf_path):
    """Generates chunks while parsing a .maf file.

    TODO: This code was extracted and modified from
    _remove_duplicate_alignments(). However the original code was not changed
    to use this new procedure.

    Yields object with keys:
        * type (string) - ('HEADER' or 'ALIGNMENT')
        * chunk (string) - multiple lines starting with score.
    """

    # Use state machine strategy.
    class READ_STATE:
        READING_HEADER = 0
        READING_ALIGNMENTS = 1

    state = READ_STATE.READING_HEADER
    current_alignment_chunk = ''

    with open(input_maf_path) as input_fh:
        for line in input_fh:
            if state == READ_STATE.READING_HEADER and line[0] != 'a':
                yield {
                    'type': 'HEADER',
                    'chunk': line
                }
            elif line[0] == 'a':
                # Permanently change state. Remaining lines come in chunks
                # representing an alignment.
                state = READ_STATE.READING_ALIGNMENTS

                # Arrived at start of alignment chunk.
                if current_alignment_chunk:
                    yield {
                        'type': 'ALIGNMENT',
                        'chunk': current_alignment_chunk
                    }
                # Start new alignment chunk.
                current_alignment_chunk = line
            else:
                # Append line to accumulating chunk.
                current_alignment_chunk += line

    # Yield the last chunk.
    yield {
        'type': 'ALIGNMENT',
        'chunk': current_alignment_chunk
    }


def build_in_memory_conservation_structure(input_maf_path):
    """Parses a .maf file and builds an in-memory structure that maps from
    positions in MG1655 to alignment blocks.
    """
    start_pos_to_data_map = {}

    chunk_generator = _generate_maf_chunks(input_maf_path)
    while True:
        try:
            chunk_data = chunk_generator.next()

            if chunk_data['type'] == 'HEADER':
                continue

            lines = chunk_data['chunk'].split('\n')

            # Get rid of blank lines:
            lines = [l for l in lines if l.strip()]

            # Grab MG1655 position.
            try:
                start_pos = _get_ecoli_position_from_maf_line(lines[1])
            except:
                print chunk_data['chunk']
                assert False

            data = {}
            for line in lines[1:]:
                parsed_line = _parse_alignment_line(line)
                data[parsed_line['strain']] = parsed_line
            assert not start_pos in start_pos_to_data_map, (
                    "Unintended overwrite.")
            start_pos_to_data_map[start_pos] = data

        except StopIteration:
            break

    return start_pos_to_data_map


class ChunkFinder(object):

    def __init__(self, start_pos_to_data_map):
        self.start_pos_to_data_map = start_pos_to_data_map
        self.sorted_start_pos_list = sorted(self.start_pos_to_data_map.keys())

        self.corresponding_end_pos_list = []
        for chunk_start in self.sorted_start_pos_list:
            chunk_data = self.start_pos_to_data_map[chunk_start]
            chunk_size = chunk_data[MG1655_MAF_STRAIN_NAME]['size']
            self.corresponding_end_pos_list.append(chunk_start + chunk_size)


    def find_chunks_for_interval(self, interval, offset=0):
        """Returns chunks that contain interval relative to MG1655 reference
        frame.

        Args:
            interval: Tuple. 0-indexed. Start and end inclusive.
                E.g. (1000,1002) is 3-base interval at position 1000.
        """
        chunks = []

        # Update interval to match offset.
        interval = tuple([p + offset for p in interval])

        # Look between these indeces. Narrowed down below.
        min_idx = 0
        max_idx = len(self.sorted_start_pos_list) - 1

        # Iteratively narrow down range of where to look.
        narrowed_down = False
        while not narrowed_down:
            assert max_idx >= min_idx, "min: %d, max: %d" % (min_idx, max_idx)
            mid_idx = (max_idx + 1 + min_idx) / 2
            if mid_idx == max_idx:
                narrowed_down = True
            chunk_start = self.sorted_start_pos_list[mid_idx]
            chunk_end = self.corresponding_end_pos_list[mid_idx]

            if chunk_end < interval[0]:
                # Not far enough.
                min_idx = mid_idx
                continue

            if interval[1] < chunk_start:
                # Too far.
                max_idx = mid_idx
                continue

            # Otherwise, we've closed in on relevant region.
            narrowed_down = True

        # Strategy: Initial brute force strategy.
        # Check each block sequentially until we find the right one.
        for start_pos_idx in range(min_idx, max_idx + 1):
            chunk_start = self.sorted_start_pos_list[start_pos_idx]
            chunk_end = self.corresponding_end_pos_list[start_pos_idx]

            # Not far enough.
            #         >    <
            # xxxxxx
            if chunk_end < interval[0]:
                continue

            # Too far.
            # >    <
            #         xxxxxx
            if interval[1] < chunk_start:
                break

            chunks.append(self.start_pos_to_data_map[chunk_start])

        # assert len(chunks)
        return chunks

    def get_synonyms_for_interval_in_chunk(self, interval, chunk_data,
            offset=0):
        """Returns dict mapping from strain name to synonym for interval
        in chunk.
        """
        synonym_data = {}

        # Update interval to match offset.
        interval = tuple([p + offset for p in interval])

        # First we figure out which indeces in the chunk we want by using
        # the MG1655 positions.
        full_mg1655_seq = chunk_data[MG1655_MAF_STRAIN_NAME]['seq']
        start_offset = interval[0] - chunk_data[MG1655_MAF_STRAIN_NAME]['start']
        interval_size = interval[1] - interval[0] + 1

        # Count how many '-' occur before offset position.
        blank_count = 0
        base_count = 0
        interval_seq = ''
        interval_base_indices = []
        for base_idx in range(len(full_mg1655_seq)):
            base = full_mg1655_seq[base_idx]
            if base == '-':
                blank_count += 1
                continue
            base_count += 1
            if base_count > start_offset:
                interval_seq += base
                interval_base_indices.append(base_idx)
            if len(interval_seq) == interval_size:
                break

        # Get other positions.
        chunk_strain_keys = chunk_data.keys()
        for strain in chunk_strain_keys:
            strain_seq = chunk_data[strain]['seq']
            strain_interval_seq = ''.join([strain_seq[base_idx]
                    for base_idx in interval_base_indices])
            synonym_data[strain] = strain_interval_seq

        return synonym_data


AGR_SET = set(['AGA', 'AGG'])

AGR_SYNONYM_SET = set(['CGA', 'CGC', 'CGG', 'CGT'])


def debug_find_subs_for_AGR(genbank_path):
    """Aggregates stats on how many AGR positions have alternates across
    the genome.
    """
    seq_record = get_genome_record(genbank_path)
    result = replace_AGR_in_genome_using_conservation(seq_record)
    print 'total', result['total_AGRs']
    print 'found', result['found']


def find_agr_intervals(feature, seq_record):
    """Returns inclusive, 0-indexed intervals.
    """
    global_frame_interval_list = []

    # Some features are composed of a JOIN of multiple sub-features. We can
    # handle this generally by leveraging interface that gives locations
    # as list.
    for loc in feature.location.parts:
        feature_seq = str(loc.extract(seq_record.seq))

        # Aggregate raw list of intervals local to to feature.
        interval_list = []
        for base_pos in range(0, len(feature_seq), 3):
            codon = feature_seq[base_pos:base_pos + 3]
            if codon in AGR_SET:
                interval_list.append((base_pos, base_pos + 2))

        # If negative-strand, update positions.
        if loc.strand == -1:
            interval_list = [(len(loc) - 1 - i[1], len(loc) - 1 - i[0])
                    for i in interval_list]

        # Update positions to global frame.
        global_frame_interval_list += [(i[0] + loc.start, i[1] + loc.start)
                for i in interval_list]

    # Assert that each interval found actually maps to AGR.
    for i in global_frame_interval_list:
        codon = str(seq_record.seq[i[0]:i[1] + 1])
        if feature.strand == -1:
            codon = reverse_complement(codon)
        assert codon in AGR_SET, (str(i), str(interval_list), codon,
                feature.qualifiers['gene'][0], feature.location.start)

    # Make 1-indexed
    # interval_list = [(i[0] + 1, i[1] + 1) for i in interval_list]

    return global_frame_interval_list


def replace_AGR_in_genome_using_conservation(seq_record,
        ignore_translation_gene_list=[], verbose=False):
    """Replaces as many AGRs in genome as possible using conservation.

    Returns dictionary of results:
        seq_record: SeqRecord with AGR replaced where possible.
    """
    # Final SeqRecord to be returned.
    updated_seq_record = SeqRecord("")

    # Conservation data.
    start_pos_to_data_map = build_in_memory_conservation_structure(
            ECOLI_CONSERVATION_MAF_UPDATED)
    chunk_finder = ChunkFinder(start_pos_to_data_map)

    # Keep track of metadata for summary.
    total_AGR_occurrences = 0
    num_replaced = 0

    # Use divider strategy.
    divider = Divider(seq_record)
    index = 0
    for div in divider.get_next_division():
    # for div in [wholediv]:
        print ('found division from', index, 'to', div,
                'total AGR', total_AGR_occurrences,
                'num replaced', num_replaced)

        # # DEBUG: Break early.
        # if index > 50000:
        #     updated_seq_record += seq_record[index:]
        #     break

        # Keep copy of current_partition in case recode fails.
        current_partition = seq_record[index:div]
        recoded_partition = copy.deepcopy(current_partition)

        # Attempt recode.
        partition_coding_features = [f for f in recoded_partition.features
                if f.type == 'CDS']
        agr_occurrences_in_partition = 0
        num_repalced_in_partition = 0
        for f in partition_coding_features:
            result = replace_AGR_in_feature_using_conservation(
                    f, recoded_partition, chunk_finder, offset=index)
            update_feature_seq_given_feature(recoded_partition, f,
                    result['updated_seq'], allow_translation_change=True)
            for annotation_f in result['codon_swap_annotations']:
                add_feature_to_seq_record(recoded_partition, annotation_f)
            agr_occurrences_in_partition += result['AGR_occurrences']
            num_repalced_in_partition += result['num_replaced']
            print f.location.start, agr_occurrences_in_partition, num_repalced_in_partition
        try:
            check_all(current_partition, recoded_partition,
                    special_cases=ignore_translation_gene_list)
            total_AGR_occurrences += agr_occurrences_in_partition
            num_replaced += num_repalced_in_partition
            updated_seq_record += recoded_partition
        except:
            updated_seq_record += current_partition

        # Update iteration index.
        index = div

    return {
        'seq_record': updated_seq_record,
        'total_AGRs': total_AGR_occurrences,
        'found': num_replaced,
    }


def replace_AGR_in_feature_using_conservation(f, seq_record, chunk_finder,
        offset=0):
    """Creates new sequence for feature with AGR codons replaced using
    conservation data. If no suitable replacement is found, AGR codon may be
    left in place.

    Returns dictionary with keys:
        * updated_seq
        * AGR_occurrences
        * num_replaced
    """
    results_dict = {}

    original_f_seq = str(seq_record.seq[f.location.start:f.location.end])
    updated_seq = copy.copy(original_f_seq)

    # New annotations indicating swaps.
    codon_swap_annotations = []

    interval_list = find_agr_intervals(f, seq_record)
    num_synonyms_found = 0
    for i in interval_list:
        # Get all conservation data chunks matching this interval.
        chunks = chunk_finder.find_chunks_for_interval(i, offset=offset)

        # If no chunks for interval, continue to next interval.
        if not chunks:
            continue

        # Else try to find synonyms for chunk.
        replacement_codon = None
        for chunk_data in chunks:
            if replacement_codon is not None:
                # No need to look at other chunks.
                break

            synonym_data = (
                    chunk_finder.get_synonyms_for_interval_in_chunk(
                            i, chunk_data, offset=offset))
            for strain, seq in synonym_data.iteritems():
                if strain == MG1655_MAF_STRAIN_NAME:
                    continue

                # Check against synonym set of correct polarity.
                if f.strand == 1:
                    # Just make sure feature has strand value.
                    pass
                elif f.strand == -1:
                    seq = reverse_complement(seq)
                else:
                    raise AssertionError("Unknown feature strand.")
                if seq in AGR_SYNONYM_SET:
                    replacement_codon = seq
                    break

        # Replace codon if found. Update count.
        if replacement_codon is not None:
            local_interval = (i[0] - f.location.start,
                    i[1] - f.location.start)
            current_codon = updated_seq[local_interval[0]:local_interval[1] + 1]
            if f.strand == -1:
                current_codon = reverse_complement(current_codon)
            assert current_codon in AGR_SET, (current_codon, i, local_interval)
            assert seq in AGR_SYNONYM_SET

            if f.strand == -1:
                seq = reverse_complement(seq)

            updated_seq = (
                    updated_seq[:local_interval[0]] +
                    seq +
                    updated_seq[local_interval[1] + 1:])

            # Use global frame of this partition for position of feature.
            start = i[0]
            new_f = SeqFeature(
                    location=FeatureLocation(start, start + 3),
                    type=CONSERVATION_CODON_SWAP_TYPE,
                    strand=f.strand,
                    id=CONSERVATION_CODON_SWAP_TYPE + '_' + str(uuid4())[:8]
            )
            codon_swap_annotations.append(new_f)

            # Update found count.
            num_synonyms_found += 1

    if f.strand == -1:
        updated_seq = reverse_complement(updated_seq)

    results_dict['updated_seq'] = updated_seq
    results_dict['codon_swap_annotations'] = codon_swap_annotations
    results_dict['AGR_occurrences'] = len(interval_list)
    results_dict['num_replaced'] = num_synonyms_found

    return results_dict


if __name__ == '__main__':
    # update_mg1655_positions(
    #         ECOLI_CONSERVATION_MAF_ORIGINAL, ECOLI_CONSERVATION_MAF_UPDATED)
    # start_pos_to_data_map = build_in_memory_conservation_structure(
    #         ECOLI_CONSERVATION_MAF_ORIGINAL)
    # # find_chunks_for_interval((103999, 104002), start_pos_to_data_map)
    # # find_chunks_for_interval((105211, 105214), start_pos_to_data_map)
    # find_chunks_for_interval((106568, 106571), start_pos_to_data_map)
    debug_find_subs_for_AGR(MG1655_GENBANK)
