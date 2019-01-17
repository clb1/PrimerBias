#!/usr/bin/env python

# Standardizing read order
# ------------------------
# 1) Place both reads in their forward orientation
# 2) Mark locations of 5p, 3p, UDTD5, and UDTD7 in each read in its forward orientation
# 3) Determine which read is the upstream read and which is the downstream read
# 4) Trim the upstream read from its 5p. Return 5p UMI
# 5) Trim the upstream read from its 3p, if it contains a 3p. Return 3p UMI.
# 6) Trim the downstream read from its 3p. Return 3p UMI.
# 7) Trim the downstream read from its 5p, if it contains a 5p. Return 5p UMI.
# 8) Reverse complement the downstream read
# 9) Error-correct the UMIs
# 10) Output: 5p UMI, upstream read, revcomp downstream read, 3p UMI
#


from Bio import AlignIO, SeqIO
from collections import defaultdict, Counter
from copy import copy
import gzip
from operator import itemgetter
import os
import re
from string import maketrans, translate
import subprocess
import sys
import pysam

import pdb

# Sequencing reads are kept or placed in their 5'->3' orientation for the strand on which they occur
class SeqRead(object):
    DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

    def __init__(self, ID, nucleotide_sequence, quality_string, quality_scores, from_reverse_strand):
        self.ID = ID
        self.original_sequencing_readID = ID
        self.seq = nucleotide_sequence
        self.qual_string = quality_string
        self.qual_scores = quality_scores

        # Start and stop positions are in the context of the 5' and 3' position of the primer/sequencing adapter molecules.
        # So if a primer is sense with the sequencing read, start < stop. But if antisense, start > stop.
        self.orientation_5p = None
        self.start_5p = None
        self.stop_5p = None
        #self.umi_5p = None
        #self.umi_5p_qual = None

        self.orientation_3p = None
        self.start_3p = None
        self.stop_3p = None
        #self.umi_3p = None
        #self.umi_3p_qual = None

        self.orientation_udtd5 = None
        self.stop_udtd5 = None
        self.start_udtd5 = None

        self.orientation_udtd7 = None
        self.stop_udtd7 = None
        self.start_udtd7 = None        

        if (from_reverse_strand):
            self.reverse_complement()

            
    def resetID(self, new_ID):
        self.ID = new_ID
        

    def getID(self):
        return self.ID


    def getSeqLen(self):
        return len(self.seq)


    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result


    def reverse_complement(self):
        self.seq = translate(self.seq, self.DNA_complement_table)[::-1]
        self.qual_string = self.qual_string[::-1]
        self.qual_scores = self.qual_scores[::-1]


    def getSequence(self, as_reverse_complement=False):
        ret_seq = self.seq
        if (as_reverse_complement):
            ret_seq = translate(self.seq, self.DNA_complement_table)[::-1]
        return ret_seq


    def getQualString(self, return_reversed=False):
        ret_string = self.qual_string
        if (return_reversed):
            ret_string = self.qual_string[::-1]
        return ret_string


    def getAdapterTrimmedSequence(self):
        reject_reason = None
        
        if (self.orientation_udtd5 != None and self.orientation_udtd7 != None): # Should only contain one adapter, at most
            reject_reason = "Read has both adapters"
            ret_seq = ""
        else:
            ret_seq = self.seq

            try:
                if (self.orientation_udtd5 != None):
                    assert (self.orientation_udtd5 == "antisense")
                    ret_seq = ret_seq[0:self.stop_udtd5]
                elif (self.orientation_udtd7 != None):
                    assert (self.orientation_udtd7 == "antisense")
                    ret_seq = ret_seq[0:self.stop_udtd7]
            except AssertionError:
                reject_reason = "Read has 'sense' adapter"
                
        return ret_seq, reject_reason

    
    def asFasta(self, as_reverse_complement=False):
        assert (len(self.seq) > 0 and len(self.qual_string)>0)
        ret_string = None
        if (as_reverse_complement):
            revcomp_seq = translate(self.seq, self.DNA_complement_table)[::-1]
            ret_string = ">%s\n%s" % (self.ID, revcomp_seq)
        else:
            ret_string = ">%s\n%s" % (self.ID, self.seq[::-1])
        return ret_string

    
    def asFastq(self, as_reverse_complement=False):
        assert (len(self.seq) > 0 and len(self.qual_string)>0)
        ret_string = None
        if (as_reverse_complement):
            revcomp_seq = translate(self.seq, self.DNA_complement_table)[::-1]
            ret_string = "@%s\n%s\n+\n%s" % (self.ID, revcomp_seq, self.qual_string[::-1])
        else:
            ret_string = "@%s\n%s\n+\n%s" % (self.ID, self.seq, self.qual_string)
        return "@%s\n%s\n+\n%s" % (self.ID, self.seq, self.qual_string)

    
    def addMatch(self, match_spec, match_start_pos, match_stop_pos, start_pos_in_reversed_coords):
        """Record the match as it occurs in the sequencing read's relative context (ie native 5'->3' orientation)."""
        assert (match_spec in ["5p","3p", "UDTD5", "UDTD7"])
        reject_reason = None

        if (match_start_pos < 0 or match_stop_pos > len(self.seq)):
            reject_reason = "Read starts in a primer"
        else:
            if (start_pos_in_reversed_coords):
                orientation = "antisense"
                match_start_pos, match_stop_pos = (len(self.seq)-match_start_pos, len(self.seq)-match_stop_pos)
            else:
                orientation = "sense"

            if (match_spec == "5p"):
                if (self.start_5p != None):
                    reject_reason = "More than two primers in a mate"
                else:
                    self.orientation_5p = orientation
                    self.start_5p = match_start_pos
                    self.stop_5p = match_stop_pos
                    #self.umi_5p = umi
                    #self.umi_5p_qual = umi_qual
            elif (match_spec == "3p"):
                if (self.start_3p != None):
                    reject_reason = "More than two primers in a mate"
                else:
                    self.orientation_3p = orientation
                    self.start_3p = match_start_pos
                    self.stop_3p = match_stop_pos
                    #self.umi_3p = umi
                    #self.umi_3p_qual = umi_qual
            elif (match_spec == "UDTD5" and self.orientation_udtd5 == None):
                self.orientation_udtd5 = orientation
                self.start_udtd5 = match_start_pos
                self.stop_udtd5 = match_stop_pos
            elif (match_spec == "UDTD7" and self.orientation_udtd7 == None):
                self.orientation_udtd7 = orientation
                self.start_udtd7 = match_start_pos
                self.stop_udtd7 = match_stop_pos

        return reject_reason

    
    def hasPrimer(self, primer_spec, orientation):
        assert (orientation in ["sense", "antisense"])
        assert (primer_spec in ["5p", "3p"])
        return ((primer_spec == "5p" and self.orientation_5p == orientation) or
                (primer_spec == "3p" and self.orientation_3p == orientation))


    def hasAPrimer(self):
        return (self.orientation_5p != None or self.orientation_3p != None)
    

    def trimToPrimers(self):
        assert (self.start_5p != None or self.start_3p != None)
        reject_reason = None

        #print >> sys.stderr, "Before trim: %s" % self.seq
        if (self.orientation_5p == "sense"):
            start = self.start_5p
            if (self.orientation_3p != None):
                if (self.orientation_3p == "antisense"):
                    stop = self.start_3p
                else:
                    reject_reason = "Both primers in same read have same sense"
            else:
                stop = len(self.seq)
        else: 
            assert (self.orientation_5p == "antisense")
            stop = self.start_5p
            if (self.orientation_3p != None):
                if (self.orientation_3p == "sense"):
                    start = self.start_5p
                else:
                    reject_reason = "Both primers in same read have same sense"
            else:
                start = 0
            
        if (reject_reason == None):
            if (stop <= start):
                reject_reason = "Antisense primer occurs before sense primer in read"
            else:
                self.seq = self.seq[start:stop]
                self.qual_string = self.qual_string[start:stop]
                self.qual_scores = self.qual_scores[start:stop]
                #print >> sys.stderr, "After  trim: %s\n" % self.seq

        return reject_reason


    def numPrimersFound(self):
        return int(self.start_5p != None) + int(self.start_3p != None)

    
def readPrimerSeqs(primers_fasta, template_name):
    all_primer_seqs = {}

    primers_5p, primers_3p = {}, {}
    for record in SeqIO.parse(primers_fasta,'fasta'):
        assert (record.id.startswith(template_name))
        if (record.id.endswith("_5p") or record.id.endswith("_3p")): # For skipping Illumina adapters
            primer_label, fwd_or_rev = record.id.rsplit('_',1)
            if (fwd_or_rev == "5p"):
                primers_5p[primer_label] = record.seq
            elif (fwd_or_rev == "3p"):
                primers_3p[primer_label] = record.seq

    for primer_label, fwd in primers_5p.items():
        fwd_rc = fwd.reverse_complement()
        rev = primers_3p[primer_label]
        rev_rc = rev.reverse_complement()
        all_primer_seqs[primer_label] = (str(fwd), str(rev), str(fwd_rc), str(rev_rc))

    return all_primer_seqs


def readTemplateSeq(templates_fasta, template_name):
    template_seq = None
    for record in SeqIO.parse(templates_fasta,'fasta'):
        if (record.id == template_name):
            assert (template_seq == None), "Assertion Error: read unexpected second template sequence."
            template_seq = str(record.seq)
    assert (template_seq != None), "Assertion Error: didn't read a template sequence."

    degenerate_region_lens = tuple(map(len, filter(lambda x: x!='', re.split("[ACGT]+", template_seq))))
    return template_seq, degenerate_region_lens

#
# upstream read:
#  o Has sense 5p primer on and then revcomp 3p on fwd strand
#
# downstream read:
#  o Has antisense 3p primer on fwd strand
def standardizeMateOrder(mate1, mate2):
    reject_reason = None
    # Standardize here to make mate1 be the mate with the 5' primer.
    mate1_has_sense_5p = mate1.hasPrimer("5p", "sense")
    mate1_has_sense_3p = mate1.hasPrimer("3p", "sense")
    mate1_has_antisense_5p = mate1.hasPrimer("5p", "antisense")
    mate1_has_antisense_3p = mate1.hasPrimer("3p", "antisense")

    mate2_has_sense_5p = mate2.hasPrimer("5p", "sense")
    mate2_has_sense_3p = mate2.hasPrimer("3p", "sense")
    mate2_has_antisense_5p = mate2.hasPrimer("5p", "antisense")
    mate2_has_antisense_3p = mate2.hasPrimer("3p", "antisense")

    # PrimerBias: First mate should have sense 5p and antisense 3p. Second mate should have only antisense 5p.

    # TODO
    # -If has revcomp 5p/revcomp 3p, then must have 3p/5p. Use for for recovery
    #if (reject_reason == None):
    #    reject_reason = assessAttemptRecovery(mate1, mate2)

    #if (mate1_has_sense_5p and mate2_has_sense_3p):
    #    if (mate1_has_antisense_3p or mate2_has_antisense_5p): # Both should have the other primer if one does
    #        if (not (mate1_has_antisense_3p and mate2_has_antisense_5p)):
    #            reject_reason = "Second primer not detected in read(s)"
    #elif (mate2_has_sense_5p and mate1_has_sense_3p):
    #    if (mate2_has_antisense_3p or mate1_has_antisense_5p): # Both should have the other primer if one does
    #        if (not (mate2_has_antisense_3p and mate1_has_antisense_5p)):
    #            reject_reason = "Second primer not detected in read(s)"
    #    mate1, mate2 = mate2, mate1
    #elif (mate1.numPrimersFound() != mate2.numPrimersFound()):
    #    reject_reason = "Mates have unequal number of primer matches"
    if (mate2_has_sense_5p and mate2_has_antisense_3p and mate1_has_antisense_5p):
        mate1, mate2 = mate2, mate1
    elif (not(mate1_has_sense_5p and mate1_has_antisense_3p and mate2_has_antisense_5p)):
        reject_reason = "Uncategorized standardizeMateOrder() failure"

    return mate1, mate2, reject_reason


def findMatchStartStop(match_spec, bam_line, reference_length):
    reject_reason = None
    found_alen_mismatch = False

    try:
        match_start_pos = bam_line.qstart
        match_stop_pos = bam_line.qend
    except SystemError:
        aligned_pairs = filter(lambda x: x[1] != None, bam_line.get_aligned_pairs())
        match_start_pos = aligned_pairs[0][0]
        match_stop_pos = aligned_pairs[-1][0]
        
    if (match_spec[-2:] in ["5p","3p"] and bam_line.alen != reference_length):
        # Fix cases where there is a 1bp on the primer's 3' end or up to a 2bp mismatch on the primer's 5' end
        aligned_pairs = filter(lambda x: x[1] != None, bam_line.get_aligned_pairs())
        if (aligned_pairs[0][1] <= 2): # 5' end
            match_start_pos -= aligned_pairs[0][1]
        elif (aligned_pairs[-1][1] == reference_length-2): # 3' end
            match_stop_pos += 1
        else:
            found_alen_mismatch = True

    elif (match_spec in ["UDTD5", "UDTD7"] and bam_line.alen != reference_length):
        # Fix cases where there is up to a 2bp mismatch on the sequencing adapter's 3' end
        aligned_pairs = filter(lambda x: x[1] != None, bam_line.get_aligned_pairs())
        if (aligned_pairs[-1][1] in [reference_length-3, reference_length-2]):
            match_stop_pos += aligned_pairs[-1][1]

    return (match_start_pos, match_stop_pos, found_alen_mismatch)


# One mate should have sense 5p and antisense 3p. The other mate should have antisense 5p.
def checkForPrimerDimer(mate1, mate2, primer_seqs):
    reject_reason = None
    aux_info = None
    
    fwd_primer, rev_primer, revcomp_fwd_primer, revcomp_rev_primer = primer_seqs

    len_fwd_primer = len(fwd_primer)
    len_rev_primer = len(rev_primer)
    
    trimmed_read1_seq, reject_reason = mate1.getAdapterTrimmedSequence()
    if (reject_reason == None):
        trimmed_read2_seq, reject_reason = mate2.getAdapterTrimmedSequence()
    
    if (reject_reason == None):
        mate1_has_sense_5p = mate1.hasPrimer("5p", "sense")
        mate1_has_sense_3p = mate1.hasPrimer("3p", "sense")
        mate1_has_antisense_5p = mate1.hasPrimer("5p", "antisense")
        mate1_has_antisense_3p = mate1.hasPrimer("3p", "antisense")

        mate2_has_sense_5p = mate2.hasPrimer("5p", "sense")
        mate2_has_sense_3p = mate2.hasPrimer("3p", "sense")
        mate2_has_antisense_5p = mate2.hasPrimer("5p", "antisense")
        mate2_has_antisense_3p = mate2.hasPrimer("3p", "antisense")

        mate1_only_5p = (mate1_has_sense_5p or mate1_has_antisense_5p) and (not (mate1_has_sense_3p or mate1_has_antisense_3p))
        mate2_only_5p = (mate2_has_sense_5p or mate2_has_antisense_5p) and (not (mate2_has_sense_3p or mate2_has_antisense_3p))

        mate1_only_3p = (not (mate1_has_sense_5p or mate1_has_antisense_5p)) and (mate1_has_sense_3p or mate1_has_antisense_3p)
        mate2_only_3p = (not (mate2_has_sense_5p or mate2_has_antisense_5p)) and (mate2_has_sense_3p or mate2_has_antisense_3p)

        mate1_has_both_improper = ((mate1_has_sense_5p and mate1_has_sense_3p) or (mate1_has_antisense_3p and mate1_has_antisense_5p))
        mate2_has_both_improper = ((mate2_has_sense_5p and mate2_has_sense_3p) or (mate2_has_antisense_3p and mate2_has_antisense_5p))

        mate1_has_both_proper = ((mate1_has_sense_5p and mate1_has_antisense_3p) or (mate1_has_sense_3p and mate1_has_antisense_5p))
        mate2_has_both_proper = ((mate2_has_sense_5p and mate2_has_antisense_3p) or (mate2_has_sense_3p and mate2_has_antisense_5p))

    return reject_reason, aux_info


def processBamLinesForTarget(target_seq_ID, bam_lines, references, reference_lengths):
    reject_reason = None
    mate1, mate2 = None, None
    any_alen_mismatch_found = False

    # Order the bam lines so that primary and better matches are first
    bam_lines = map(lambda x: (x, int(x.seq!=None), x.get_tag("AS"), int(not(x.is_supplementary or x.is_secondary))), bam_lines)
    bam_lines = sorted(bam_lines, key=itemgetter(1,2,3), reverse=True)
    
    for bam_line in map(lambda y: y[0], bam_lines):
        found_alen_mismatch = False
        match_spec = references[bam_line.tid]
        reference_length = reference_lengths[bam_line.tid]
        
        if (bam_line.is_unmapped):
            reject_reason = "At least one mate unmapped"
            break

        # TODO: check UDTD alen correction
        match_start_pos, match_stop_pos, found_alen_mismatch = findMatchStartStop(match_spec, bam_line, reference_length)
        
        start_pos_in_reversed_coords = bam_line.is_reverse
        any_alen_mismatch_found = any_alen_mismatch_found or found_alen_mismatch

        # Match is to *_5p, *_3p, UDTD5, or UDTD7
        if (match_spec not in ["UDTD5", "UDTD7"]):
            assert (target_seq_ID == references[bam_line.tid][0:-3])
            match_spec = references[bam_line.tid][-2:]  

        if (bam_line.qname.endswith("_1")):
            if (mate1 == None):
                mate1 = SeqRead(bam_line.qname[0:-2], bam_line.seq, bam_line.qual, bam_line.query_qualities, bam_line.is_reverse)
            reject_reason = mate1.addMatch(match_spec, match_start_pos, match_stop_pos, start_pos_in_reversed_coords)
        elif (bam_line.qname.endswith("_2")):
            if (mate2 == None):
                mate2 = SeqRead(bam_line.qname[0:-2], bam_line.seq, bam_line.qual, bam_line.query_qualities, bam_line.is_reverse)
            reject_reason = mate2.addMatch(match_spec, match_start_pos, match_stop_pos, start_pos_in_reversed_coords)

        if (reject_reason != None):
            break
        
    if (reject_reason == None):
        assert(mate1 != None and mate2 != None)
        if (not (mate1.hasAPrimer() and mate2.hasAPrimer())):
            reject_reason = "Mate w/o primer match"
            mate1, mate2 = None, None
        
    if (reject_reason == None and any_alen_mismatch_found):
        reject_reason = "Incomplete primer match"
        
    return (reject_reason, mate1, mate2)


def trimAndStandardizeMates(target_seq_ID, bam_lines, references, reference_lengths, primer_seqs):
    """1) Trim reads upstream of primer sequence(s), relative to primer sequences(s)
       2) Relabel reads <umi_5p>_<umi_3p>
    """
    aux_info = None
    reject_reason, mate1, mate2 = processBamLinesForTarget(target_seq_ID, bam_lines, references, reference_lengths)

    if (reject_reason == None):
        mate1, mate2, reject_reason = standardizeMateOrder(mate1, mate2)

    #
    # TODO: Replace functionality
    #
    if (reject_reason == None):
        reject_reason = mate1.trimToPrimers()
        if (reject_reason == None):
            mate2.trimToPrimers()
            # TODO: not sure if this is relevant now
            #if (reject_reason == None):
            #    reject_reason = performMatesQC(mate1, mate2, primer_seqs, umi_len)

    #    umis_5p = mate1.getUMI("5p") + mate2.getUMI("5p")
    #    umis_5p_qual = mate1.getUMIQual("5p") + mate2.getUMIQual("5p")
    #    umis_3p = mate1.getUMI("3p") + mate2.getUMI("3p")
    #    umis_3p_qual = mate1.getUMIQual("3p") + mate2.getUMIQual("3p")
    #
    #    reject_reason, umi_based_ID = formUMIbasedPairID(target_seq_ID, umis_5p, umis_5p_qual, umis_3p, umis_3p_qual, umi_len)
    #
    #    if (reject_reason == None):
    #        mate1.resetID(umi_based_ID)
    #        mate2.resetID(umi_based_ID)
    #
        
    return (reject_reason, aux_info, mate1, mate2)

        
def performMatesQC(mate1, mate2, primer_seqs, umi_len):
    """To be called after order of the mates has been standardized and the mates' sequences
    have been trimmed to their primer starts"""
    reject_reason = None
    
    fwd_primer, rev_primer, revcomp_fwd_primer, revcomp_rev_primer = primer_seqs

    read1_seq = mate1.getSequence()
    read2_seq = mate2.getSequence()

    len_fwd_primer = len(fwd_primer)
    len_rev_primer = len(revcomp_rev_primer)
    sum_primer_lens = len_fwd_primer + len_rev_primer

    if (len(read1_seq) < sum_primer_lens or len(read2_seq) < sum_primer_lens):
        reject_reason = "QC fail: trimmed mate too short"
    else:
        read1_starts_w_fwd_primer = len([i for i in xrange(len_fwd_primer) if read1_seq[i] != fwd_primer[i]]) <= 4
        read2_starts_w_rev_primer = len([i for i in xrange(len_rev_primer) if read2_seq[i] != rev_primer[i]]) <= 4
        
        if (not read1_starts_w_fwd_primer or not read2_starts_w_rev_primer):
            reject_reason = "QC fail: primer doesn't initiate read"
            
    return reject_reason


def mergeMatesAndExtractSequences(mate1, mate2, template_seq):
    reject_reason, UMI, primer_region = None, None, None

    mate1_seq = mate1.getSequence()
    mate1_qual = mate1.getQualString()
    mate2_seq = mate2.getSequence(True)
    mate2_qual = mate2.getQualString(True)

    three_seqs = ">template\n%s\n>mate1\n%s\n>mate2\n%s" % (template_seq, mate1_seq, mate2_seq)
    child = subprocess.Popen("muscle -quiet", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    child.stdin.write(three_seqs)
    child.stdin.close()
    align = AlignIO.read(child.stdout, "fasta")
    align.sort() # Places the alignments in the order mate1, mate2, template

    # Merging heuristics:
    #
    # Cases t=='N'
    # ------------
    # m1 == m2        -> m1
    # m1 != m2        -> higher qval
    # m1 == '-'       -> m2
    # m2 == '-'       -> m1
    # m1 & m2 == '-'  -> skip to make shorter UMI. Should be rare

    # Case t=='-' within a degenerate region
    # -----------
    # m1 or m2 == '-' -> skip, erroneous
    # m1 == m2 != '-' -> add m1  to make longer UMI. Should be rare
    # m1 != m2        -> higher qval, (or skip?). Should be rare
    # m1 & m2 == '-'  -> shouldn't happen

    # Logic for within degenerate region, based on template sequence nucleotides
    # prev_nuc  curr_nuc   in_degenerate
    #   nuc       nuc         no
    #   nuc       N           yes
    #   nuc       -           no
    #   N         nuc         no
    #   N         N           yes
    #   N         -           yes if next template nuc is 'N'
    #   -         nuc         no
    #   -         N           yes
    #   -         -           state stays the same

    S = ''
    m1_ind, m2_ind, t_ind = -1, -1, -1
    prev_t = 'A' # dummy nucleotide to start
    # Add a nucleotide to S at every degenerate region positions
    for c in xrange(len(align[0])):
        m1,m2,t = list(align[:,c])
        if (m1 != '-'): m1_ind += 1  # These three are indices into the nucleotide sequences, not into MSA positions. That would be 'c'.
        if (m2 != '-'): m2_ind += 1
        if (t != '-'):   t_ind += 1                

        if (t != '-'):
            in_degenerate_region = t == 'N'
        else:
            if (prev_t == 'N'):
                in_degenerate_region = t_ind+1<len(template_seq) and template_seq[t_ind+1]=='N'
            elif (prev_t != '-'):
                in_degenerate_region = False
            
        if (t == 'N'):
            if (m1 == m2 == '-'):
                S += ' '
            elif (m1 == '-'):
                S += m2
            elif (m2 == '-'):
                S += m1
            elif (m1 == m2):
                S += m1
            else:
                S += m1 if (mate1_qual[m1_ind] >= mate2_qual[m2_ind]) else m2
        elif (t == '-'):
            if (m1 == '-' or m2 == '-'):
                assert (m1 != '-' or m2 != '-')
                S += ' '
            elif (m1 == m2):
                S += m1 if (in_degenerate_region) else '.'
            elif (m1 != m2):
                S += ' '
        else:
            S += '.'
            
        prev_t = t

    msa = "### %s ###\n" % mate1.getID()
    msa += "%s\tmate1\n" % str(align[0].seq)
    msa += "%s\tmate2\n" % str(align[1].seq)
    msa += "%s\ttemplate\n" % str(align[2].seq)
    msa += "%s\tdegenerate region seqs\n\n" % S

    if (False):
        print >> sys.stderr, msa

    degenerate_region_seqs = tuple(filter(lambda y: y!='', map(lambda x: x.replace(' ',''), re.split("\.+", S))))

    if (len(degenerate_region_seqs) < 2):
        reject_reason = "Missing degenerate region"
        degenerate_region_seqs = (None, None)
    elif (len(degenerate_region_seqs) > 2):
        reject_reason = "Too many degenerate regions"
        degenerate_region_seqs = (None, None)

    return (reject_reason, msa, degenerate_region_seqs)


def processReadPairAlignments(curr_readID_bam_lines, references, reference_lengths, reject_reason_target_counts, all_primer_seqs):
    reject_reason, aux_info = None, None
    mate1, mate2 = None, None
    target_seq_ID = None
    
    primer_targets = set(map(lambda y: y[0:-3], filter(lambda x: x.endswith("_5p") or x.endswith("_3p"), curr_readID_bam_lines.keys())))

    if (len(primer_targets) > 1):
        # Multiple different target region primers matched to this one read pair
        key = tuple(sorted(primer_targets))
        reject_reason = "Primers for different targets in same read pair"
    elif (len(primer_targets)==0):
        alignment_targets = set(curr_readID_bam_lines.keys())
        if (alignment_targets.issubset(set(["UDTD5", "UDTD7"]))):
            reject_reason = "Only adapter sequences"
        else:
            bam_lines = []
            for l in curr_readID_bam_lines.values():
                bam_lines.extend(l)
            assert(all(map(lambda x: x.is_unmapped, bam_lines)))
            reject_reason = "Read pair unmapped"
    else:
        target_seq_ID = list(primer_targets)[0]
        bam_lines = []
        for l in curr_readID_bam_lines.values():
            bam_lines.extend(l)
        primer_seqs = all_primer_seqs[target_seq_ID]
        reject_reason, aux_info, mate1, mate2 = trimAndStandardizeMates(target_seq_ID, bam_lines, references, reference_lengths, primer_seqs)

        if (reject_reason != None):
            reject_reason_target_counts[reject_reason][target_seq_ID] += 1
            
    return (reject_reason, aux_info, target_seq_ID, mate1, mate2)


def modifyReadPairsUsingAlignments(sam_stream, output_tsv, output_aux, all_primer_seqs, template_seq, degenerate_region_lens):
    reject_reason_statistics = defaultdict(int)
    reject_reason_target_counts = defaultdict(lambda: defaultdict(int))

    ip_bam = pysam.AlignmentFile(sam_stream, 'r')
    op = gzip.open(output_tsv, 'wb')
    op.write("# Sequences are from the priming strand, written 5'->3'\n")
    op.write("# Source_Read_Pair\tInternal_UMI\tPriming_Sequence\n")
    op_aux = gzip.open(output_aux, 'wb')
    
    references = ip_bam.references
    reference_lengths = ip_bam.lengths
    curr_readID_bam_lines = defaultdict(list)
    num_read_pairs_read, num_read_pairs_written = 0,0
    num_read_pairs_w_no_primer_matches, num_read_pairs_w_incomplete_primer_matches = 0,0

    bam_line = ip_bam.next()

    curr_readID = bam_line.qname[0:-2]
    curr_target_seq_ID = ip_bam.references[bam_line.tid]
    curr_readID_bam_lines[curr_target_seq_ID].append(bam_line)

    try:
        while (True):
            next_bam_line = ip_bam.next()
            next_readID = next_bam_line.qname[0:-2]
            curr_target_seq_ID = references[next_bam_line.tid]

            if (next_readID != curr_readID):
                reject_reason, aux_info, target_seq_ID, mate1, mate2 = processReadPairAlignments(curr_readID_bam_lines, references, reference_lengths,
                                                                                                 reject_reason_target_counts, all_primer_seqs)
                num_read_pairs_read += 1
                
                if (reject_reason == None):
                    reject_reason, msa, (UMI, primer_region) = mergeMatesAndExtractSequences(mate1, mate2, template_seq)

                if (reject_reason == None and ((len(UMI), len(primer_region)) != degenerate_region_lens)):
                    reject_reason = "Incorrect length of degenerate region(s)"
                    #print >> sys.stderr, (UMI, primer_region)
                    #pdb.set_trace()

                if (reject_reason == None):
                    op.write("%s\t%s\t%s\n" % (curr_readID, UMI, primer_region))
                    num_read_pairs_written += 1

                reject_reason_statistics[reject_reason] += 1

                if (reject_reason != None and aux_info != None):
                    op_aux.write("%s\t%s\t%s\n" % (reject_reason, target_seq_ID, aux_info))
                    
                curr_readID_bam_lines.clear()
                curr_readID_bam_lines[curr_target_seq_ID].append(next_bam_line)
                curr_readID = next_readID

            else:
                curr_readID_bam_lines[curr_target_seq_ID].append(next_bam_line)

    except StopIteration:
        reject_reason, aux_info, target_seq_ID, mate1, mate2 = processReadPairAlignments(curr_readID_bam_lines, references, reference_lengths,
                                                                                         reject_reason_target_counts, all_primer_seqs)
        num_read_pairs_read += 1
        reject_reason_statistics[reject_reason] += 1
        
        if (reject_reason == None):
            reject_reason, msa, (UMI, primer_region) = mergeMatesAndExtractSequences(mate1, mate2, template_seq)

        if (reject_reason == None):
            op.write("%s\t%s\t%s\n" % (target_seq_ID, UMI, primer_region))
            num_read_pairs_written += 1

        if (reject_reason != None and aux_info != None):
            op_aux.write("%s\t%s\t%s\n" % (reject_reason, target_seq_ID, aux_info))

    ip_bam.close()
    op.close()
    op_aux.close()

    return (reject_reason_statistics, num_read_pairs_read, num_read_pairs_written, reject_reason_target_counts)


def writeLogfile(num_read_pairs_read, num_read_pairs_written, reject_reason_statistics, reject_reason_target_counts):
    consistency_count = 0
    
    op = open(logfile, 'w')
    op.write("Total number of read pairs as input: %d\n" % num_read_pairs_read)

    perc_written = 100.0 * float(num_read_pairs_written)/float(num_read_pairs_read)
    op.write("Number of good read pairs written: %d (%4.2f%%)\n" % (num_read_pairs_written, perc_written))
    consistency_count += num_read_pairs_written
    
    rejections_to_list_out = []
    reason_count_tups = reject_reason_statistics.items()
    reason_count_tups = sorted(reason_count_tups, key=itemgetter(1), reverse=True)
    op.write("\nRejection Statistics:\n")
    for reason, count in reason_count_tups:
        if (reason != None):
            perc = 100.0 * float(count)/float(num_read_pairs_read)
            op.write("%s: %d (%4.2f%%)\n" % (reason, count, perc))
            consistency_count += count
            if (perc > 1.0):
                rejections_to_list_out.append( (perc, reason) )

    if (consistency_count != num_read_pairs_read):
        op.write("\nWARNING: Cannot account for all of the read pairs read. Consistency count = %d.\n" % consistency_count)

    op.close()


if (__name__ == "__main__"):
    template_name, primers_fasta, templates_fasta, sam_stream, output_tsv, output_aux, logfile = sys.argv[1:]

    all_primer_seqs = readPrimerSeqs(primers_fasta, template_name)
    template_seq, degenerate_region_lens = readTemplateSeq(templates_fasta, template_name)
    
    (reject_reason_statistics, num_read_pairs_read, num_read_pairs_written, reject_reason_target_counts) = \
      modifyReadPairsUsingAlignments(sam_stream, output_tsv, output_aux, all_primer_seqs, template_seq, degenerate_region_lens)

    writeLogfile(num_read_pairs_read, num_read_pairs_written, reject_reason_statistics, reject_reason_target_counts)

    sys.exit(0)
    
