#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Project: trichCode
# Created by dong on  12:35 2014/12/20
# Version $

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import subprocess
import shutil
import shelve
import MySQLdb as mdb

VERSION = '1.2.0'
DESCRIPTION = 'Run trichCode genome annotation pipeline'

KEYS = dict()

# region UTILITIES

def _get_cmd_dirname(tool_name):
    return subprocess.check_output('dirname `which %s`' % tool_name, shell=True)


def _get_dir_parent(dirct):
    return '/' + '/'.join(dirct.split('/')[1:len(dirct.split('/')) - 1])


def is_fasta(in_f):
    try:
        with open(in_f) as f:
            if not f.readline().startswith('>'):
                return False
    except Exception, e:
        print >> sys.stderr, e
    return True 


def _init_shelve_db():
    global KEYS

    print os.path.join(KEYS['DIR_project'], 'stat.db')
    s = shelve.open(os.path.join(KEYS['DIR_project'], 'stat.db'))
    s.close()


def _set_mark_stat(key, stat):
    global KEYS

    s = shelve.open(os.path.join(KEYS['DIR_project'], 'stat.db'))
    try:
        s[key] = stat
    finally:
        s.close()


def _get_mark_stat(key):
    global KEYS

    s = shelve.open(os.path.join(KEYS['DIR_project'], 'stat.db'))
    try:
        if s.has_key(key):
            stat = s[key]
        else:
            stat = None
    finally:
        s.close()
    return stat


def _get_contig_name(seq):
    for line in open(seq):
        line = line.strip()
        if line.startswith('>'):
            return line[1:5]


def check_depend_files(list_files):
    global KEYS
    
    for f in list_files:
        # check file exist
        if KEYS.has_key(f) and os.path.isfile(KEYS[f]):
            continue
        else:
            # This need to be optimised in future
            print >> sys.stderr, "required file is missing: " + KEYS[f]
            return False
    return True

# endregion of Utilities ------------------------


# region PASA -----------------------------------

def run_pasa_init_assm():
    """
    before run pasa, initialize MySQL database
    :return:
    """
    global KEYS
    print '[Initializing MySQL database]'
    try:
        con = mdb.connect('localhost', KEYS['mysql_user'], KEYS['mysql_passwd'])
        cur = con.cursor()
        cur.execute('DROP DATABASE %s' % KEYS['name_project'])
    except mdb.Error, e:
        print 'Error %d: %s' % (e.args[0], e.args[1])
    finally:
        if con:
            con.close()

    # copy pasa.alignAssembly.config.cfg to project folder and replace the db name
    cmd_3 = 'cd %s/transcripts; cp %s/configs/pasa.alignAssembly.Template.txt alignAssembly.config.cfg;' \
            'sed -i \'s/<__MYSQLDB__>/%s/g\' alignAssembly.config.cfg; ' \
            'sed -i \'s/<__MIN_PERCENT_ALIGNED__>/%s/g\' alignAssembly.config.cfg; ' \
            'sed -i \'s/<__MIN_AVG_PER_ID__>/%s/g\' alignAssembly.config.cfg;' % (
                KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['name_project'], KEYS['min_percent_aligned'], KEYS['min_avg_per_id'])
    
    print '[start to run: PASA init]\n' + cmd_3
    subprocess.check_call(cmd_3, shell=True)


def run_pasa_init_makeup():
    """
    before run pasa, initialize MySQL database
    :return:
    """
    global KEYS

    # copy pasa.annotationCompare.config.cfg to project folder and replace the db name
    cmd = 'cd %s/transcripts; cp %s/configs/pasa.annotationCompare.Template.txt annotationCompare.config.cfg;' \
            ' sed -i \'s/<__MYSQLDB__>/%s/g\' annotationCompare.config.cfg' % (
                KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['name_project'])
    print '[start to run: PASA makeup]\n' + cmd
    subprocess.check_call(cmd, shell=True)

def run_pasa_align_assem():
    '''
    alignment assembly
    :produce: [name_project].pasa_assemblies.gff3
    '''
    global KEYS

    stat = _get_mark_stat('pasa_assembly')
    if stat is not None and os.path.isfile(stat):
        KEYS['gff_pasa_assemblies'] = stat
        return
    
    run_pasa_init_assm()

    name_transcript = KEYS['DIR_project'] + '/transcripts/transcript.fas'
    shutil.copyfile(KEYS['seq_transcript'], name_transcript)
    print KEYS['dir_pasa']
    cmd_1 = 'cd %s/transcripts; %s/seqclean/seqclean/seqclean %s > %s/seqclean.log 2>' \
            ' %s/seqclean.stderr' % (KEYS['DIR_project'], KEYS['dir_pasa'], name_transcript, KEYS['dir_log'], KEYS['dir_log'])
    print '[start to run: PASA seqclean]\n' + cmd_1
    subprocess.check_call(cmd_1, shell=True)

    # unmasked genome is used here for obtaining more robust transcript alignment
    cmd_2 = 'cd %s/transcripts; %s/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config.cfg -C -R -I %s' \
            ' -g %s -t transcript.fas.clean -T -u transcript.fas --ALIGNERS blat,gmap --MAX_INTRON_LENGTH %s' \
            ' --CPU %s > %s/pasa_align_assem.log 2> %s/pasa_align_assem.stderr' % (
                KEYS['DIR_project'], KEYS['dir_pasa'], KEYS['intron_max'], KEYS['seq_genome'], KEYS['intron_max'],
                KEYS['cores'], KEYS['dir_log'], KEYS['dir_log'])
    print '[start to run: Launch_PASA_pipeline]\n' + cmd_2
    subprocess.check_call(cmd_2, shell=True)
    KEYS['gff_pasa_assemblies'] = KEYS['DIR_project'] + '/transcripts/%s.pasa_assemblies.gff3' % KEYS['name_project']
    _set_mark_stat('pasa_assembly', KEYS['gff_pasa_assemblies'])

    print KEYS['gff_pasa_assemblies']

def run_pasa_extract_ORFs():
    '''
    Extraction of ORFs from PASA assemblies (auto-annotation and/or reference ORFs for training gene predictors)
    :produce: [name_project].assemblies.fasta.transdecoder.genome.gff3
    '''
    global KEYS
    
    stat = _get_mark_stat('pasa_extract')
    if stat is not None and os.path.isfile(stat):
        KEYS['gff_pasa_training'] = stat
        return
    
    cmd = 'cd %s/transcripts; %s/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta %s --pasa_transcripts_gff3 %s ' \
          '> %s/pasa_extract_ORFs.log 2> %s/pasa_extract_ORFs.stderr' % (
              KEYS['DIR_project'], KEYS['dir_pasa'], KEYS['name_project'] + '.assemblies.fasta',
              KEYS['name_project'] + '.pasa_assemblies.gff3', KEYS['dir_log'], KEYS['dir_log'])
    print '[start to run: transcode]\n' + cmd
    subprocess.check_call(cmd, shell=True)
    KEYS['gff_pasa_training'] = KEYS['DIR_project'] + '/transcripts/%s.assemblies.fasta.transdecoder.genome.gff3' % \
                                                        KEYS['name_project']
    _set_mark_stat('pasa_extract', KEYS['gff_pasa_training'])

def run_pasa_update_gene_model():
    '''
    Update, to incorporate the PASA alignment evidence, correcting exon boundaries, adding UTRs, and
    models for alternative splicing based on the PASA alignment assemblies generated above.
    :return:
    '''
    global KEYS

    stat = _get_mark_stat('gff_pasa_update')
    if stat is not None and os.path.isfile(stat):
        KEYS['gff_pasa_update'] = stat
        return
    
    run_pasa_init_makeup()
    
    cmd_1 = 'cd %s/transcripts; %s/scripts/Launch_PASA_pipeline.pl -c %s/transcripts/annotationCompare.config.cfg' \
            ' -g %s -t %s/transcripts/transcript.fas.clean -A -L --annots_gff3 %s/predicts/evm/evm.gff3' \
            ' > %s/pasa_update_gene_model.log 2> %s/pasa_update_gene_model.stderr' % (
                KEYS['DIR_project'], KEYS['dir_pasa'], KEYS['DIR_project'], KEYS['seq_genome'],
                KEYS['DIR_project'], KEYS['DIR_project'], KEYS['dir_log'], KEYS['dir_log'])
    print '[start to run: PASA update]\n' + cmd_1
    subprocess.check_call(cmd_1, shell=True)

    # find the updated gff3, copy it to products folder
    cmd_2 = 'cd %s/transcripts; cp `ls %s.*post*gff3 -1t | head -n 1` %s/products/pasa_update.gff3 ' % (
        KEYS['DIR_project'], KEYS['name_project'], KEYS['DIR_project'])
    print '[start to run: PASA upate rename]\n' + cmd_2
    subprocess.check_call(cmd_2, shell=True)
    KEYS['gff_pasa_update'] = KEYS['DIR_project'] + '/products/pasa_update.gff3'

    cmd_0 = 'cd %s/products/; python %s/scripts/statistics.py pasa_update.gff3 pasa > %s/stat_pasa_update.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_0
    subprocess.check_call(cmd_0, shell=True)
    
    _set_mark_stat('gff_pasa_update', KEYS['gff_pasa_update'])
# endregion


# region PREPARE TRAINING SET
def run_transposonPSI():
    global KEYS

    stat = _get_mark_stat('tpsi')
    if stat is not None and os.path.isfile(stat):
        return

    cmd = 'cd %s; transposonPSI.pl %s nuc > %s/transposonPSI.log 2> %s/transposonPSI.stderr' % (
        KEYS['DIR_project'] + '/' + 'repeats/transposon', KEYS['seq_genome'], KEYS['dir_log'],
        KEYS['dir_log'])
    print '[start to run: transposonPSI]\n' + cmd
    # transposonPSI.pl target_test_genome_seq.fasta nuc
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call('pwd', shell=True)

    subprocess.check_call('cd %s;mv *TPSI.allHits.chains.bestPerLocus.gff3 %s' % (
        KEYS['DIR_project'] + '/' + 'repeats/transposon', KEYS['gff_transposon']), shell=True)

    _set_mark_stat('tpsi', KEYS['gff_transposon'])


def run_GeneMarkET():
    """
    :produce file:
    """
    global KEYS

    stat = _get_mark_stat('gmet')
    if stat is not None and os.path.isfile(stat):
        return

    print stat
    print 'should not be here'
    exit(255)

    et_score = 0
    if KEYS['gff_intron_gmet'] != '':
        cmd = 'cd %s;gmes_petap.pl --sequence %s -ET %s --fungus --et_score %s ' \
              '--cores %s > %s/GeneMark-ET.log 2> %s/GeneMark-ET.stderr' % (
                  KEYS['DIR_project'] + '/' + 'predicts/gm', KEYS['seq_genome_masked'],
                  KEYS['gff_intron_gmet'], KEYS['score_intron'], KEYS['cores'], KEYS['dir_log'],
                  KEYS['dir_log'])
    else:
        # cmd = 'cd %s;gmes_petap.pl %s > %s/GeneMark-ET.log 2> %s/GeneMark-ET.stderr' \
        # % (KEYS['DIR_project'] + '/' + 'predicts/gm', KEYS['seq_genome_masked'],
        # KEYS['dir_log'], KEYS['dir_log'])
        cmd = 'cd %s;gmes_petap.pl --sequence %s -ES --fungus ' \
              '--cores %s > %s/GeneMark-ET.log 2> %s/GeneMark-ET.stderr' % (
                  KEYS['DIR_project'] + '/' + 'predicts/gm', KEYS['seq_genome_masked'],
                  KEYS['cores'], KEYS['dir_log'], KEYS['dir_log'])
    print '[start to run: GeneMarkET]\n' + cmd
    subprocess.check_call(cmd, shell=True)

    _set_mark_stat('gmet', KEYS['gtf_gm'])
    return


def run_exonerate():
    global KEYS

    stat = _get_mark_stat('exonerate')
    if stat is not None and os.path.isfile(stat):
        return

    cmd = 'cd %s; exonerate --maxintron %s --minintron %s --model p2g --showvulgar yes ' \
          '--showalignment no --showquerygff no --showtargetgff yes --percent %s ' \
          '--ryo "AveragePercentIdentity: %%pi\\n\" %s %s >> %s 2>> %s/exonerate.stderr' % (
              KEYS['DIR_project'], KEYS['intron_max'], KEYS['intron_min'], KEYS['exon_similarity'],
              KEYS['seq_homology'], KEYS['seq_genome'], KEYS['gff_homology'], KEYS['dir_log'])
    print '[start to run: Exonerate]\n' + cmd
    subprocess.check_call(cmd, shell=True)

    _set_mark_stat('exonerate', KEYS['gff_homology'])


# find match
def run_find_match(filter_rule):
    global KEYS
    if KEYS.has_key('gff_pasa_assemblies') and KEYS['gff_pasa_assemblies'] != '':
        cmd_1 = 'cd %s; python %s/scripts/find_match.py --transposon %s --transcriptalg %s --transcripttype pasa --exonerate %s ' \
                '--genemark-et %s --out matched_set.gtf 2>> %s/find_matched.stderr' % (
                    KEYS['DIR_project'] + '/' + 'training', KEYS['DIR_trichoCODE'],
                    KEYS['gff_transposon'], KEYS['gff_pasa_training'], KEYS['gff_homology'],
                    KEYS['gtf_gm'], KEYS['dir_log'])

        print '[start to run: find_match.py with transcript data]\n' + cmd_1
        subprocess.check_call(cmd_1, shell=True)
    else:
        cmd_1 = 'cd %s; python %s/scripts/find_match.py --transposon %s --exonerate %s ' \
                '--genemark-et %s --out matched_set.gtf 2>> %s/find_matched.stderr' % (
                    KEYS['DIR_project'] + '/' + 'training', KEYS['DIR_trichoCODE'],
                    KEYS['gff_transposon'], KEYS['gff_homology'], KEYS['gtf_gm'], KEYS['dir_log'])

        print '[start to run: find_match.py without transcript data]\n' + cmd_1
        subprocess.check_call(cmd_1, shell=True)

    cmd_2 = 'cd %s/training; python %s/scripts/statistics.py matched_set.gtf match > %s/stat_matched_set.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_2
    subprocess.check_call(cmd_2, shell=True)

    cmd_3 = 'cd %s; python %s/scripts/prepare_training.py --matched_set %s --filter_rule %s --training_set %s >>' \
            '%s/prepare_training.log 2>> %s/prepare_training.stderr' % (
                KEYS['DIR_project'] + '/' + 'training', KEYS['DIR_trichoCODE'],
                KEYS['gtf_matched_set'], filter_rule, KEYS['gtf_training_set'], KEYS['dir_log'],
                KEYS['dir_log'])
    print '[start to run: prepare_training.py]\n' + cmd_3
    subprocess.check_call(cmd_3, shell=True)

    cmd4 = 'cd %s/training/; %s/scripts/gtf2gff3.pl training_set.gtf > training_set.gff3; ' % (
               KEYS['DIR_project'], KEYS['DIR_trichoCODE'])
    print '[EVM transfomating Genemark-ET output to gff3]:'
    subprocess.check_call(cmd4, shell=True)

    cmd_5 = 'cd %s/training; python %s/scripts/statistics.py training_set.gff3 gff >> %s/stat_matched_set.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_5
    subprocess.check_call(cmd_5, shell=True)
    

# endregion prepare training set  ---------------------------


# region TRAINING

def _aug_training():
    global KEYS
    """
    Training augustus
    :return: None
    """
    cmd_0 = 'cd %s; %s/scripts/gff2gbSmallDNA.pl %s %s 5000 training.gbk --overlap >> %s/gff2gbSmallDNA.log ' \
            ' 2>> %s/gff2gbSmallDNA.stderr' % (
        KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['dir_augustus'],
        KEYS['gtf_training_set'], KEYS['seq_genome'], KEYS['dir_log'], KEYS['dir_log'])
    subprocess.check_call(cmd_0, shell=True)
    subprocess.check_call(
        'cd %s; grep -c LOCUS training.gbk* > %s/AUGUSTUS_trainingset.stat.log' % (
            KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['dir_log']), shell=True)

    cmd_2 = 'rm %s -rf; %s/scripts/new_species.pl --species=%s' % (
        KEYS['dir_augustus'] + '/' + 'config/species/' + KEYS['name_project'], KEYS['dir_augustus'], \
        KEYS['name_project'])
    subprocess.check_call(cmd_2, shell=True)

    cmd_3 = 'cd %s; etraining --species=%s training.gbk >> %s/etraining.log 2>> %s/etraining.stderr' % (
        KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['name_project'], KEYS['dir_log'],
        KEYS['dir_log'])
    subprocess.check_call(cmd_3, shell=True)


def _aug_est_hints():
    global KEYS
    """
    :return: hints.est.gff
    """
    cmd_1 = 'cd %s; blat -noHead -maxIntron=%s %s %s est.psl > %s/blat.log 2> %s/blat.stderr' % (
        KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['intron_max'], KEYS['seq_genome'],
        KEYS['seq_transcript'], KEYS['dir_log'], KEYS['dir_log'])
    subprocess.check_call(cmd_1, shell=True)

    cmd_2 = 'cd %s; cat est.psl | %s/scripts/filterPSL.pl --best --minCover=80 > est.f.psl ' \
            '2> %s/filterPSL.stderr' % (KEYS['DIR_project'] + '/' + 'predicts/augustus',
         KEYS['dir_augustus'], KEYS['dir_log'])
    subprocess.check_call(cmd_2, shell=True)

    cmd_3 = 'cd %s; %s/scripts/blat2hints.pl --nomult --in=est.f.psl --out=hints.est.gff 2>' \
            ' %s/blat2hints.stderr' % (KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['dir_augustus'],
                KEYS['dir_log'])
    subprocess.check_call(cmd_3, shell=True)

    return '%s/hints.est.gff' % (KEYS['DIR_project'] + '/' + 'predicts/augustus')


def _aug_intron_hints():
    if os.path.isfile(KEYS['gff_intron_aug']):
        return KEYS['gff_intron_aug']


def _aug_repeat_hints():
    """
    :return: hints.est.gff
    """
    global KEYS
    cmd = 'cat %s  | tail -n +3 | perl -ne \'chomp; s/^\s+//; @t = split(/\s+/); ' \
          'print $t[4].\"\t\".\"repmask\tnonexonpart\t\".$t[5].\"\t\".$t[6].\"\t0\t.\t.\tsrc=RM\n\";\' ' \
          '| sort -n -k 1,1 > repeats.gff 2> %s/prepare_repeat_hints.stderr; sed \'1d\' repeats.gff > temp; mv temp repeats.gff' % (
              KEYS['repeat_out'], KEYS['dir_log'])
    subprocess.check_call(cmd)
    return '%s/repeats.gff' % (KEYS['DIR_project'] + '/' + 'predicts/augustus')


def _aug_homology_hints():
    """
    :return: hints.homology.gff
    """
    global KEYS
    cmd = 'cd %s; %s/scripts/exonerate2hints.pl --in=%s --source=P --out=exonerate.hints' % (
        KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['dir_augustus'], KEYS['gff_homology'])
    subprocess.check_call(cmd, shell=True)
    return '%s/exonerate.hints' % (KEYS['DIR_project'] + '/' + 'predicts/augustus')


def _aug_prepare_hints():
    global KEYS

    print 'preparing hints for AUGUSTUS ...  '
    hints_total = list()
    if os.path.isfile(KEYS['gff_homology']):
        print 'hint from homology proteins'
        hints_hom = _aug_homology_hints()
        hints_total.append(hints_hom)
    if os.path.isfile(KEYS['seq_transcript']):
        print 'hint from transcript sequences'
        hints_est = _aug_est_hints()
        hints_total.append(hints_est)
    if os.path.isfile(KEYS['gff_intron_aug']):
        print 'hint from intron'
        hints_intron = _aug_intron_hints()
        hints_total.append(hints_intron)

    for hint in hints_total:
        cmd = 'cd %s; rm %s -f; cat %s >> %s' % (
            KEYS['DIR_project'] + '/' + 'predicts/augustus', 'hints.gff', hint, 'hints.gff')
        subprocess.check_call(cmd, shell=True)


# endregion


# region PREDICTION

def run_AUGUSTUS():
    global KEYS

    stat = _get_mark_stat('augustus')
    if stat is not None and os.path.isfile(stat):
        KEYS['gff_augustus'] = stat
        return

    cmd = 'cd %s; augustus --species=%s %s --extrinsicCfgFile=%s --hintsfile=hints.gff > augustus.hints.gff' \
          % (KEYS['DIR_project'] + '/' + 'predicts/augustus', KEYS['name_project'],
        KEYS['seq_genome_masked'], KEYS['DIR_trichoCODE'] + '/configs/extrinsic.cfg')
    print '[start to run: Augustus]\n' + cmd
    subprocess.check_call(cmd, shell=True)

    KEYS['gff_augustus'] = KEYS['DIR_project'] + '/' + 'predicts/augustus/augustus.hints.gff'
    _set_mark_stat('augustus', KEYS['gff_augustus'])


def run_SNAP():
    """
    :output file: snap_<name_project>:
    """
    global KEYS

    cmd = 'cd %s; %s/scripts/gtf2gff3.pl %s > training_set.gff3' % (
        KEYS['DIR_project'] + '/' + 'predicts/snap', KEYS['DIR_trichoCODE'],
        KEYS['gtf_training_set'])
    subprocess.check_call(cmd, shell=True)

    stat = _get_mark_stat('snap')
    if stat is not None and os.path.isfile(stat):
        KEYS['zff_snap'] = stat
        return

    cmd_1 = 'cd %s; sh %s/scripts/run_SNAP.sh %s/scripts/maker2/maker2zff ' \
            ' %s %s > %s/run_SNAP.log 2> %s/run_SNAP.stderr' % (
                KEYS['DIR_project'] + '/' + 'predicts/snap', KEYS['DIR_trichoCODE'],
                KEYS['DIR_trichoCODE'], KEYS['seq_genome'], KEYS['name_project'], KEYS['dir_log'],
                KEYS['dir_log'])
    try:
        print '[start to run: SNAP]\n' + cmd_1
        subprocess.check_call(cmd_1, shell=True)
    except subprocess.CalledProcessError, e:
        print e
        exit(255)

    KEYS['zff_snap'] = KEYS['DIR_project'] + '/' + (
        'predicts/snap/snap_%s.zff' % KEYS['name_project'])
    _set_mark_stat('snap', KEYS['zff_snap'])


# endregion


# region COMBINATION OF MULTIPLE EVIDENCES

def evm_pre_predictor():
    global KEYS

    cmd = 'cd %s/predicts/evm/; %s/EvmUtils/misc/SNAP_output_to_gff3.pl %s/predicts/snap/snap_%s.zff %s > snap.gff3;' \
          '%s/scripts/delete_mark.py \# snap.gff3 snap.gff3-2;mv snap.gff3-2 snap.gff3;' \
          '%s/EvmUtils/gff3_gene_prediction_file_validator.pl snap.gff3 > %s/evm_pre_SNAP.log 2> %s/evm_pre_SNAP.stderr' % (
              KEYS['DIR_project'], KEYS['dir_evm'], KEYS['DIR_project'], KEYS['name_project'],
              KEYS['seq_genome'], KEYS['DIR_trichoCODE'], KEYS['dir_evm'], KEYS['dir_log'],
              KEYS['dir_log'])
    print '[EVM transfomating SNAP output to gff3]:\n' + cmd
    subprocess.check_call(cmd, shell=True)

    cmd_0 = 'cd %s/predicts/evm/; python %s/scripts/statistics.py snap.gff3 gff > %s/stat_snap.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_0
    subprocess.check_call(cmd_0, shell=True)
    
    cmd1 = 'cd %s/predicts/evm/; grep -P \"AUGUSTUS\\t(CDS|exon)\" %s/predicts/augustus/augustus.hints.gff > augustus.gtf;' \
           '%s/scripts/gtf2gff3.pl augustus.gtf > augustus.gff3 2> %s/log/transformat_augustus.log;rm augustus.gtf;' \
           '%s/EvmUtils/gff3_gene_prediction_file_validator.pl augustus.gff3' % (
               KEYS['DIR_project'], KEYS['DIR_project'], KEYS['DIR_trichoCODE'],
               KEYS['DIR_project'], KEYS['dir_evm'])
    print '[EVM transfomating AUGUSTUS output to gff3]:'
    subprocess.check_call(cmd1, shell=True)

    cmd_0 = 'cd %s/predicts/evm/; python %s/scripts/statistics.py augustus.gff3 gff > %s/stat_augustus.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_0
    subprocess.check_call(cmd_0, shell=True)
    
    cmd2 = 'cd %s/predicts/evm/; %s/scripts/gtf2gff3.pl %s > genemark_hmm.gff3; ' \
           '%s/EvmUtils/gff3_gene_prediction_file_validator.pl genemark_hmm.gff3' % (
               KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['gtf_gm'], KEYS['dir_evm'])
    print '[EVM transfomating Genemark-ET output to gff3]:'
    subprocess.check_call(cmd2, shell=True)

    cmd_0 = 'cd %s/predicts/evm/; python %s/scripts/statistics.py genemark_hmm.gff3 gff > %s/stat_genemark.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_0
    subprocess.check_call(cmd_0, shell=True)
    
    cmd3 = 'cd %s/predicts/evm/; %s/scripts/delete_mark.py vulgar %s exonerate-2;mv exonerate-2 exonerate;' \
           '%s/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl exonerate > exonerate.gff3;rm exonerate;' \
           '%s/EvmUtils/gff3_gene_prediction_file_validator.pl exonerate.gff3' % (
               KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['gff_homology'], KEYS['dir_evm'],
               KEYS['dir_evm'])
    print '[EVM transfomating Exonerate output to gff3]:'
    subprocess.check_call(cmd3, shell=True)


def evm_pre_pasa():
    global KEYS
    cmd = 'cd %s/predicts/evm/; %s/EvmUtils/gff3_gene_prediction_file_validator.pl %s/transcripts/%s.pasa_assemblies.gff3;' \
          'cp %s/transcripts/%s.pasa_assemblies.gff3 pasa.gff3' % (
              KEYS['DIR_project'], KEYS['dir_evm'], KEYS['DIR_project'], KEYS['name_project'],
              KEYS['DIR_project'], KEYS['name_project'])
    subprocess.check_call(cmd, shell=True)


def evm_pre_weights():
    global KEYS
    cmd = 'cd %s/predicts/evm/; cp %s/configs/weights .; sed -i \'s/x/%s/g\' weights' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['name_project'])
    subprocess.check_call(cmd, shell=True)


def evm_run():
    global KEYS

    stat = _get_mark_stat('evm')
    if stat is not None and os.path.isfile(stat):
        KEYS['gff3_evm'] = stat
        return

    cmd_0 = 'cd %s; rm %s* here* \##gff-ver* -rf;rm evm.* -rf' % (
        KEYS['DIR_project'] + '/predicts/evm', _get_contig_name(KEYS['seq_genome']))
    subprocess.check_call(cmd_0, shell=True)

    cmd_1 = 'cd %s; cat snap.gff3 augustus.gff3 genemark_hmm.gff3 > gene_predictions.gff3;' % (
        KEYS['DIR_project'] + '/predicts/evm')
    subprocess.check_call(cmd_1, shell=True)

    if KEYS.has_key('gff_pasa_assemblies') and KEYS['gff_pasa_assemblies'] != '':
        cmd_2 = 'cd %s; %s/EvmUtils/partition_EVM_inputs.pl --genome %s --gene_predictions gene_predictions.gff3 ' \
                ' --protein_alignments exonerate.gff3 --transcript_alignments pasa.gff3 --segmentSize 100000 --overlapSize '\
                ' 10000 --partition_listing partitions_list.out > %s/evm_partitions.log 2> %s/evm_partitions.stderr' \
                % (KEYS['DIR_project'] + '/predicts/evm', KEYS['dir_evm'], KEYS['seq_genome'],
                    KEYS['dir_log'], KEYS['dir_log'])
        subprocess.check_call(cmd_2, shell=True)

        cmd_3 = 'cd %s; %s/EvmUtils/write_EVM_commands.pl --genome %s --weights %s/predicts/evm/weights ' \
                ' --gene_predictions gene_predictions.gff3 --protein_alignments exonerate.gff3 --transcript_alignments pasa.gff3 ' \
                ' --output_file_name evm.out --partitions partitions_list.out > commands.list 2> %s/write_EVM_commands.stderr' \
                % (KEYS['DIR_project'] + '/predicts/evm', KEYS['dir_evm'], KEYS['seq_genome'],
                    KEYS['DIR_project'], KEYS['dir_log'])
        subprocess.check_call(cmd_3, shell=True)
    else:
        cmd_2 = 'cd %s; %s/EvmUtils/partition_EVM_inputs.pl --genome %s --gene_predictions gene_predictions.gff3 ' \
                ' --protein_alignments exonerate.gff3 --segmentSize 100000 --overlapSize '\
                ' 10000 --partition_listing partitions_list.out > %s/evm_partitions.log 2> %s/evm_partitions.stderr' \
                % (KEYS['DIR_project'] + '/predicts/evm', KEYS['dir_evm'], KEYS['seq_genome'],
                    KEYS['dir_log'], KEYS['dir_log'])
        subprocess.check_call(cmd_2, shell=True)

        cmd_3 = 'cd %s; %s/EvmUtils/write_EVM_commands.pl --genome %s --weights %s/predicts/evm/weights ' \
                ' --gene_predictions gene_predictions.gff3 --protein_alignments exonerate.gff3 ' \
                ' --output_file_name evm.out --partitions partitions_list.out > commands.list 2> %s/write_EVM_commands.stderr' \
                % (KEYS['DIR_project'] + '/predicts/evm', KEYS['dir_evm'], KEYS['seq_genome'],
                    KEYS['DIR_project'], KEYS['dir_log'])
        subprocess.check_call(cmd_3, shell=True)

    cmd_4 = 'cd %s; %s/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log ;' \
            '%s/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out; ' \
            '%s/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome %s ' \
            '> %s/evm_running.log 2> %s/evm_running.stderr' % (
                KEYS['DIR_project'] + '/predicts/evm', KEYS['dir_evm'], KEYS['dir_evm'],
                KEYS['dir_evm'], KEYS['seq_genome'], KEYS['dir_log'], KEYS['dir_log'])
    print '[start to run: EVM]\n' + cmd_4
    subprocess.check_call(cmd_4, shell=True)

    # produce evm.gff3
    contig_name = _get_contig_name(KEYS['seq_genome'])
    cmd_5 = 'cd %s; echo %s* > list; %s/scripts/run_cat.py list; cp %s/predicts/evm/evm.gff3 %s/products/evm.gff3' % (
        KEYS['DIR_project'] + '/predicts/evm', contig_name, KEYS['DIR_trichoCODE'], KEYS['DIR_project'], KEYS['DIR_project'])
    subprocess.check_call(cmd_5, shell=True)

    cmd_0 = 'cd %s/predicts/evm/; python %s/scripts/statistics.py evm.gff3 gff > %s/stat_evm.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_0
    subprocess.check_call(cmd_0, shell=True)
    
    KEYS['gff3_evm'] = KEYS['DIR_project'] + '/predicts/evm/evm.gff3'
    _set_mark_stat('evm', KEYS['gff3_evm'])

def run_evm_to_fasta():
    global KEYS

    print '[start to run: produce fasta files on genes, cDNA, CDS and proteins]\n'
    if KEYS.has_key('gff_pasa_update') and check_depend_files(['gff_pasa_update']):
        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s prot > %s 2> %s/VI_to_proteins.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff_pasa_update'], KEYS['seq_genome'], KEYS['name_project'] + '_pasa_update_aa.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        print KEYS['gff_pasa_update']
        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s CDS > %s 2> %s/VI_to_CDS.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff_pasa_update'], KEYS['seq_genome'], KEYS['name_project'] + '_pasa_update_cds.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s gene > %s 2> %s/VI_to_gene.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff_pasa_update'], KEYS['seq_genome'], KEYS['name_project'] + '_pasa_update_gene.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s cDNA > %s 2> %s/VI_to_cDNA.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff_pasa_update'], KEYS['seq_genome'], KEYS['name_project'] + '_pasa_update_cDNA.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)
        
        return

    if KEYS.has_key('gff3_evm_filter') and check_depend_files(['gff3_evm_filter']):
        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s prot > %s 2> %s/to_proteins.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff3_evm_filter'], KEYS['seq_genome'], KEYS['name_project'] + '_evm_aa.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s CDS > %s 2> %s/to_proteins.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff3_evm_filter'], KEYS['seq_genome'], KEYS['name_project'] + '_evm_cds.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s gene > %s 2> %s/to_proteins.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff3_evm_filter'], KEYS['seq_genome'], KEYS['name_project'] + '_evm_gene.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        cmd = 'cd %s/products/; %s/EvmUtils/gff3_file_to_proteins.pl %s %s cDNA > %s 2> %s/to_proteins.stderr' % (
            KEYS['DIR_project'], KEYS['dir_evm'], KEYS['gff3_evm_filter'], KEYS['seq_genome'], KEYS['name_project'] + '_evm_cDNA.fas', KEYS['dir_log'])
        subprocess.check_call(cmd, shell=True)

        return

    print 'check file gff_pasa_update or gff3_evm_filter are not existed!!!'


# endregion of evm --------------------------------


# region of filter -----------------------------

def filter_with_rules():
    global KEYS
    
    stat = _get_mark_stat('evm_filter')
    if stat is not None and os.path.isfile(stat):
        KEYS['gff3_evm_filter'] = stat
   
    cmd = 'cd %s/products; %s/scripts/filter_with_rules.py %s %s %s > evm_filter.gff3 2> %s/filter.stderr' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['gff3_evm'], KEYS['filter_intron_max'],
        KEYS['filter_protein_min'], KEYS['dir_log'])
    KEYS['gff3_evm_filter'] = KEYS['DIR_project'] + '/products/evm_filter.gff3'
    subprocess.check_call(cmd, shell=True)
    
    cmd_0 = 'cd %s/products/; python %s/scripts/statistics.py evm_filter.gff3 gff > %s/stat_evm_filter.txt' % (
        KEYS['DIR_project'], KEYS['DIR_trichoCODE'], KEYS['dir_rpts'])
    print '[start to run: statistics.py]\n' + cmd_0
    subprocess.check_call(cmd_0, shell=True)
    
    KEYS['gff3_evm_filter'] = KEYS['DIR_project'] + '/products/evm_filter.gff3'
    _set_mark_stat('evm_filter', KEYS['gff_augustus'])

# endregion of filter -----------------------------


# region of makeup -----------------------------

def makeup():
    # UTRs, alternatively spliced models
    if check_depend_files(['seq_transcript']):
        _set_mark_stat('gff_pasa_update', '')
        run_pasa_update_gene_model()

        # statistics of updated information
          
    # translate into proteins sequences, gene, CDS
    run_evm_to_fasta()

# endregion of makeup -----------------------------


def run_trichoCODE(manual=False, debug=False):
    run_pre_training(manual, debug)

def run_pre_training(manual=False, debug=False):
    if check_depend_files(['seq_transcript']):
        # when in debug mode, keep all the stat
        if debug == False:
            _set_mark_stat('pasa_assembly', '')
            _set_mark_stat('pasa_extract', '')

        run_pasa_align_assem()
        run_pasa_extract_ORFs()

    list_depend_files = ('seq_genome', 'seq_genome_masked', 'seq_homology')
    OK = check_depend_files(list_depend_files)
    if OK:
        if debug == False:
            _set_mark_stat('tpsi', '')
            _set_mark_stat('gmet', '')
            _set_mark_stat('exonerate', '')
        run_transposonPSI()
        run_GeneMarkET()
        run_exonerate()
    else:
        print >> sys.stderr, 'Please check the Required files in the config.cfg is right'
        exit(1)
    
    run_find_match('loose')
    if not manual:
        run_training(manual, debug)

def run_training(manual=False, debug=False):
    _aug_training()
    _aug_prepare_hints()

    if not manual:
        run_prediction(manual, debug)

def run_prediction(manual=False, debug=False):
    if debug == False:
        _set_mark_stat('augustus', '')
        _set_mark_stat('snap', '')
    run_AUGUSTUS()
    run_SNAP()
    
    if not manual:
        run_combination(manual, debug)

def run_combination(manual=False, debug=False):
    evm_pre_predictor()
    if check_depend_files(['seq_transcript', ]):
        evm_pre_pasa()
    evm_pre_weights()
    if debug == False:
        _set_mark_stat('evm', '')
    evm_run()

    if not manual:
        run_makeup()

def run_makeup():
    filter_with_rules()
   
    makeup()

def mkdir_init():
    global KEYS

    subprocess.check_call('mkdir -p %s' % KEYS['DIR_project'], shell=True)
    os.chdir(KEYS['DIR_project'])
    subprocess.check_call('sh %s' % os.path.join(KEYS['DIR_trichoCODE'], 'scripts/make_dir.sh'), shell=True)
    file_structure = subprocess.check_output(
        "ls -R | grep \":$\" | sed -e \'s/:$//\' -e \'s/[^-][^\/]*\//--/g\' -e \'s/^/   /\' -e \'s/-/|/\'", shell=True)
    print('Current project direcotry is: ' + KEYS['DIR_project'])
    print('Now we have direcotries: \n' + file_structure)

def config_loader(config_file):
    global KEYS

    _config = SafeConfigParser()
    _config.read(config_file)

    for section in _config.sections():
        print 'Section:', section
        print 'Options:', _config.options(section)
        for name, value in _config.items(section):
            # parser_config.set()
            print '   %s = %s' % (name, value)
        print

    seq_hom = _config.get("Project", 'homology_protein')
    if not is_fasta(seq_hom):
        raise ValueError('Sequence file ' + seq_hom + ' does not appear to be a fasta file')

    dir_cur = os.getcwd()
    KEYS['DIR_project'] = os.path.join(dir_cur, _config.get('Project', 'project_name'))
    KEYS['DIR_trichoCODE'] = os.path.dirname(os.path.abspath(__file__))

    KEYS['dir_augustus'] = _get_dir_parent(_get_cmd_dirname('augustus'))
    KEYS['dir_snap'] = _get_cmd_dirname('snap')
    KEYS['dir_pasa'] = _get_dir_parent(_get_cmd_dirname('pasa'))
    KEYS['dir_evm'] = _get_cmd_dirname('evidence_modeler.pl').strip()
    KEYS['dir_log'] = os.path.join(KEYS['DIR_project'], 'log')
    KEYS['dir_rpts'] = os.path.join(KEYS['DIR_project'], 'reports')

    KEYS['cores'] = _config.get('Project', 'cores')
    KEYS['name_project'] = _config.get('Project', 'project_name')
    KEYS['seq_genome'] = _config.get('Project', 'genome')
    KEYS['seq_genome_masked'] = _config.get('Project', 'masked_genome')

    genome = _config.get("Project", 'genome')
    if not is_fasta(genome):
        raise ValueError('Genome file' + genome + ' does not appear to be a fasta file')
    genome_masked = _config.get("Project", 'masked_genome')
    if not is_fasta(genome_masked):
        raise ValueError(
            'Sequence file' + genome_masked + ' does not appear to be a fasta file')
    
    # Exonerate
    KEYS['seq_homology'] = _config.get('Project', 'homology_protein')
    KEYS['intron_min'] = _config.get('Project', 'min_intron_length')
    KEYS['intron_max'] = _config.get('Project', 'max_intron_length')
    KEYS['exon_similarity'] = _config.get('Project', 'similarity')
    seq_hom = _config.get("Project", 'homology_protein')
    if not is_fasta(seq_hom):
        raise ValueError('Sequence file' + seq_hom + ' does not appear to be a fasta file')

    # RNA-Seq
    KEYS['min_percent_aligned'] = _config.get('RNA-Seq', 'min_percent_aligned')
    KEYS['min_avg_per_id'] = _config.get('RNA-Seq', 'min_avg_per_id')
    KEYS['seq_transcript'] = _config.get('RNA-Seq', 'transcripts_assembled_fas')
    KEYS['gff_intron_gmet'] = _config.get('RNA-Seq', 'intron_gff_gmet')
    KEYS['gff_intron_aug'] = _config.get('RNA-Seq', 'intron_gff_aug')
    KEYS['score_intron'] = _config.get('RNA-Seq', 'intron_score')

    # PASA
    KEYS['mysql_user'] = _config.get('Project', 'mysql_user_name')
    KEYS['mysql_passwd'] = _config.get('Project', 'mysql_user_passwd')
    if KEYS['seq_transcript'] != '':
        if KEYS['mysql_user'] == '' or KEYS['mysql_passwd'] == '':
            raise ValueError('mysql user name or password cannot be blank')

    # filter rules
    KEYS['filter_intron_max'] = _config.get('filter', 'filter_intron_max')
    KEYS['filter_protein_min'] = _config.get('filter', 'filter_protein_min')

    # internal keys
    KEYS['gff_transposon'] = os.path.join(KEYS['DIR_project'],
                                          'repeats/transposon/transposon.gff3')
    KEYS['gff_homology'] = os.path.join(KEYS['DIR_project'], 'exonerate_hom.gff3')
    KEYS['gtf_gm'] = os.path.join(KEYS['DIR_project'], 'predicts/gm/genemark.gtf')
    KEYS['gtf_matched_set'] = os.path.join(KEYS['DIR_project'], 'training/matched_set.gtf')
    KEYS['gtf_training_set'] = os.path.join(KEYS['DIR_project'], 'training/training_set.gtf')
    KEYS['hints'] = ''

    mkdir_init()
    
    # shelve db initialize
    _init_shelve_db()


def controller(stage, manual=False, debug=False):
    print debug
    if stage == 'whole' or stage == 'pre_training':
        # PREPARATION FOR TRAINING
        # seq_genome, seq_hom, seq_genome_mask, intron_gff, seq_transcript
        run_trichoCODE(manual, debug)
    elif stage == 'training':
        run_pre_training(manual, debug)
    elif stage == 'prediction':
        run_prediction(manual, debug)
    elif stage == 'combination':
        run_combination(manual, debug)
    elif stage == 'makeup':
        run_makeup()
    else:
        print "Please indicate right stage to run. trichoCODE.py -h"
        sys.exit(255)


if __name__ == '__main__':
    DEVNULL = open(os.devnull, 'w')
    programs = ['gmes_petap.pl', 'snap', 'augustus', 'exonerate', 'transposonPSI.pl',
                'evidence_modeler.pl', 'pasa', 'blat']
    for prog in programs:
        if subprocess.call('which %s' % prog, shell=True, stdout=DEVNULL, stderr=DEVNULL) != 0:
            print >> sys.stderr, 'Please install the %s and then add the installing directory containing the %s to your system PATH' % (
                prog, prog)
            DEVNULL.close()
            sys.exit(255)
    DEVNULL.close()

    # parse args
    parser_arg = argparse.ArgumentParser(description=DESCRIPTION)
    parser_arg.add_argument('--config_file', '-c', dest='config_file',
                            help='Path to configuration file')
    # whole, pre_training, training, prediction, combination, makeup, filter
    parser_arg.add_argument('--stage', '-s', dest='stage', default='whole',
                            help='Choose which stage to start, default is whole, still you can choose pre_training, '\
                            ' training, prediction, combination, makeup')
    parser_arg.add_argument('--manual', '-m', action='store_true',
                            help='mark the manual mode')
    parser_arg.add_argument('--debug', '-d', action='store_true', help='set debug mode')
    parser_arg.add_argument('--version', action='version', version='%(prog)s' + VERSION)

    args = parser_arg.parse_args()

    if not len(sys.argv) > 1:
        print parser_arg.parse_args(['-h'])
        sys.exit(255)

    if not os.path.isfile(args.config_file):
        print "config file: %s is not exist, please check."
        sys.exit(255)

    config_loader(args.config_file)
    controller(args.stage, args.manual, args.debug)

