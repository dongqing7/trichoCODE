#!/usr/bin/env python

from __future__ import division
import sys
from optparse import OptionParser
from Bio.Blast import NCBIXML

class QueryResult(object):
    def __init__(self, id):
        self.id = id
        self.query_def = ""
        self.query_len = None
        self.hit_id = ""
        self.hit_def = ""
        self.hit_len = None
        self.alg_len = None
        self.query_start = None
        self.query_end = None
        self.hit_start = None
        self.hit_end = None
        self.identities = None
        self.bit_score = None
        self.e_value = None

def get_positions_in_query_def(query_def):
    """
    Positions: from_34_to_32432
    """
    # find the index of "Positoins"
    positions_index = query_def.find("Positions") + 9
    # extract the positions info
    positions = "from" + query_def[positions_index:].replace(' ', '_')
    return positions

def get_hit_pid(hit_def):
    splits = hit_def.split('|')
    hit_pid = splits[1]
    return hit_pid

def get_gi_from_hit_def(query_result):
    hit_info = "NO GI!"

    # gi|154684519|ref|YP_001419680.1| chromosomal replication initiation protein [Bacillus amyloliquefaciens FZB42]
    if query_result.hit_id.find("gi|") > -1:
        hit_info = query_result.hit_id
    if query_result.hit_def.find("gi|") > -1:
        hit_info = query_result.hit_def

    if hit_info.find("gi") == -1:
        return hit_info
    # print("hit_info: " + hit_info)
    gi = hit_info.split('|')[1]

    return gi

def get_aminoacids_annotation(hit_def):
    # <Hit_def>DNA polymerase III subunit beta [Paenibacillus polymyxa SC2]</Hit_def>
    index = hit_def.find(' [')
    annotation = ""
    if index > -1:
        annotation = hit_def[0:index]
    else:
        annotation = hit_def
    return annotation

def parse_blast_result(infile):
    """
    for blast2.2.24
    """
    dict_query_results = dict()
    list_results_order = list()

    result_handler = open(infile)

    query_count = 0
    for query_reader in NCBIXML.parse(result_handler):
        query_count += 1
        query_result = QueryResult(id=query_reader.query)
        # print(query_def)

        query_result.query_len = query_reader.query_length

        hits = query_reader.alignments

        hits_count = len(hits)  # q_reader.num_hits
        if hits_count < 1:
            query_count -= 1
            query_result.query_def = "NO Hits!"
            dict_query_results[query_result.id] = query_result
            list_results_order.append(query_result.id)
            continue

        # just process the first hit
        hit = hits[0]

        query_result.hit_len = hit.length
        query_result.hit_id = hit.hit_id
        query_result.hit_def = hit.hit_def

        hsps = hit.hsps

        hsps_count = len(hsps)

        # just process the first hsp
        hsp = hsps[0]

        query_result.query_def = "hits counts:" + str(hits_count) + "|" + "hsps counts:" + str(hsps_count)
        query_result.identities = hsp.identities
        query_result.alg_len = hsp.align_length
        query_result.query_start = hsp.query_start
        query_result.query_end = hsp.query_end
        query_result.hit_start = hsp.sbjct_start
        query_result.hit_end = hsp.sbjct_end
        query_result.bit_score = hsp.bits
        query_result.e_value = hsp.expect

        dict_query_results[query_result.id] = query_result
        list_results_order.append(query_result.id)

    return list_results_order, dict_query_results

def process_pfam_blast_results(list_results_order, dict_query_results):
    list_results = list()

    for query_id  in list_results_order:
        query_result = dict_query_results[query_id]

        if query_result.query_def == "NO Hits!":
            line = "%s\t%s\t%s\t%s\t%s\t%s" % (query_result.id, query_result.query_def, "-", "-", "-", "-")
            # list_results.append(line)
            continue
        # <Hit_def>pfam13514, AAA_27, AAA domain.  This domain is found in a number of double-strand DNA break proteins. This domain contains a P-loop motif.</Hit_def>
        pfam_info = query_result.hit_def

        line = "%s\t%s\t%s\t%s" % (query_result.id, query_result.query_len, query_result.query_def, pfam_info)
        # print(line)
        list_results.append(line)

    return list_results

def write_results_to_file(list_results, outfile):
    outfile_handler = open(outfile, 'w')

    for result in list_results:
        outfile_handler.write(result + '\n')

    print("write file " + outfile + " done!")

if __name__ == "__main__":
#    usage = "usage: %prog [option] file_to_parse, outfile or outfile_with_fit and not file"
#
#    parser = OptionParser(usage=usage)
#    parser.add_option("-i", "--infile", dest="infile", metavar="FILE")
#    parser.add_option("-c", "--cog", dest="cog", default=False, help="parse the cog results")
#    parser.add_option("-b", "--blast", dest="blast", default=False, help="parse the blast results")
#    parser.add_option("-o", "--out", dest="out_file", help="write ", help="genome file with fasta format", metavar="FILE")
#
#    (options, args) = parser.parse_args()
#    if len(args) != 3:
#        parser.error("Incorrect number of arguments")
#
#    infile = options.infile
#    if options.cog

    if len(sys.argv) != 2:
        print "usage: %s pfam.xml" % sys.argv[0]
        sys.exit(255)

    infile = sys.argv[1]
    outfile = infile.split('.')[0] + ".csv"

    parse_results = parse_blast_result(infile)
#    parse_results = parse_blast_result_all(infile)
    list_results_order = parse_results[0]
    dict_query_results = parse_results[1]

    # process for rpsblast result of PFAM
    list_results = process_pfam_blast_results(list_results_order, dict_query_results)
    write_results_to_file(list_results, outfile)

