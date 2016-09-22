#!/usr/bin/python
"""Convert various flavours of gff or gtf to gff3 that is acceptable to GBrowse2

Author: Ian Reid
May 20, 2009
"""

# Copyright (c) 2013, The Developers
# All rights reserved.

import sys, os.path

def usage():
    print >> sys.stderr,  'Usage: python $s infile.gff' % sys.argv[0]
    sys.exit(1)
    
def do_fixgff(infile,  outfile):
    outfile.write('###\n')
    currName = ''
    tr_id =0
    exonCount = 1
    transcript = {}
    names = {}
    parts = []
    for line in infile :
        fields = line.split('\t')
        if len(fields) < 9 :
            continue
        type = fields[2]
        #print type
        attribute = fields[8]
        attr = attribute.split(';')
        attrs = {}
        for nameval in attr:
            if len(nameval.strip()) < 3:
                continue
            try:
                (name,  value) = tuple(nameval.split())
        #    print name,  value
                if value.strip('"'):
                    attrs[name] = value.strip('"')
            except ValueError, e:
                print >> sys.stderr, e
        if 'proteinId' in attrs:
            line_name = attrs['proteinId']
            if 'name' in attrs and attrs['name'] not in names:
                names[attrs['name']] = line_name
        if 'transcriptId' in attrs:
            line_name = attrs['transcriptId']
            if 'name' in attrs and attrs['name'] not in names:
                names[attrs['name']] = line_name
        elif 'gene_id' in attrs:
            line_name = attrs['gene_id']
            if 'name' in attrs and attrs['name'] not in names:
                names[attrs['name']] = line_name
        elif 'name' in attrs and attrs['name'] in names:
            line_name = names[attrs['name']]
        else:
            line_name = None
            print >> sys.stderr, 'Unrecognized feature name'
        if line_name != currName :
            if 'seqId' in transcript :
                ## output a gene record
                gene_id = 'gene' + tr_id
                transcript['attr'] = 'ID=%s;Name=%s' % (gene_id,  currName)
                if 'alias' in attrs:
                    transcript['attr'] += ';Alias=%s' % attrs['alias']
                # transcript['type'] = 'transcript'
                transcript['type'] = 'gene'
                record = '\t'.join([transcript['seqId'],  transcript['source'],  transcript['type'],  str(transcript['start']),  str(transcript['end']),  transcript['score'],  transcript['strand'], transcript['phase'],  transcript['attr']]) + '\n'
                outfile.write(record)
                
                ## output a mRNA record
                mRNA_id = 'mRNA' + tr_id
                transcript['attr'] = 'ID=%s;Name=%s;Parent=%s' % (mRNA_id,  currName, gene_id)
                transcript['type'] = 'mRNA'
                record = '\t'.join([transcript['seqId'],  transcript['source'],  transcript['type'],  str(transcript['start']),  str(transcript['end']),  transcript['score'],  transcript['strand'], transcript['phase'],  transcript['attr']]) + '\n'
                outfile.write(record)
    #        for (key,  val) in transcript.iteritems():
    #            print key,  '=',  val
                for record in parts:
                    outfile.write(record)
                outfile.write('###\n')
            parts = []
#            print 'new Name: ',  line_name
            currName = line_name
            transcript.clear()
            exonCount = 1
            CDSCount = 1
            tr_id = currName
            transcript['type'] = 'mRNA'
            transcript['phase'] = '.'
            transcript['source'] = 'fixgff'
            if 'transcriptId' in attrs :
                tr_id = attrs['transcriptId']
            elif 'transcript_id' in attrs :
                tr_id = attrs['transcript_id']
        if   'start' not in transcript or int(fields[3] < transcript['start']):
            transcript['start'] = int(fields[3])
        if  'end' not in transcript or int(fields[4]) > transcript['end']:
            transcript['end'] = int(fields[4])
        if   'seqId' not in transcript :
            transcript['seqId'] = fields[0]
        if   'score' not in transcript :
            transcript['score'] = fields[5]
        if   'strand' not in transcript :
            transcript['strand'] = fields[6]
        attrs['Parent'] = 'mRNA' + tr_id
        if 'name' in attrs:
            attrs['alias'] = attrs['name']
            del attrs["name"]
        attrs['Name'] = currName
        if fields[2] == 'exon' :
            attrs['ID'] = attrs['Parent'] + '.%s%s' % (fields[2] ,  exonCount)
            exonCount += 1
        elif fields[2] == 'CDS' :
            attrs['ID'] = attrs['Parent'] + '.%s%s' % (fields[2] ,  CDSCount)
            CDSCount += 1
        else:
            attrs['ID'] = attrs['Parent'] + '.%s' % (fields[2],)
        attrStr=''
        
        for (key,  val) in attrs.iteritems() :
            attrStr+= '%s=%s;' % (key,  val)
        fields[8]    = attrStr.rstrip(';')
        record = '\t'.join(fields) + '\n'
        parts.append(record)

    if 'seqId' in transcript :
            ## output a gene record
            gene_id = 'gene' + tr_id
            transcript['attr'] = 'ID=%s;Name=%s' % (gene_id,  currName)
            transcript['type'] = 'gene'
            record = '\t'.join([transcript['seqId'],  transcript['source'],  transcript['type'],  str(transcript['start']),  str(transcript['end']),  transcript['score'],  transcript['strand'], transcript['phase'],  transcript['attr']]) + '\n'
            outfile.write(record)
            
            ## output a mRNA record
            mRNA_id = 'mRNA' + tr_id
            transcript['attr'] = 'ID=%s;Name=%s;Parent=%s' % (mRNA_id,  currName, gene_id)
            transcript['type'] = 'mRNA'
            record = '\t'.join([transcript['seqId'],  transcript['source'],  transcript['type'],  str(transcript['start']),  str(transcript['end']),  transcript['score'],  transcript['strand'], transcript['phase'],  transcript['attr']]) + '\n'
            outfile.write(record)
            for record in parts:
                outfile.write(record)
            outfile.write('###\n')

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] == '-h':
        usage()
    infilename =  sys.argv[1]
    infile = open(infilename,'r')
    basename, extension = os.path.splitext(infilename)
    outfilename = basename + '.gff3'
    outfile = open(outfilename,'w')

    do_fixgff(infile,  outfile)
    
    outfile.flush()
    outfile.close()
    infile.close()
