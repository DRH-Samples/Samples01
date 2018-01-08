#!/usr/bin/env python

'''
Input a BED file. Sort it, convert it to bigBed, copy the bigBed to a specified tracks folder,
and update tracks configuration using specified names, as follow:
bigBed filename = <trackName>.bb
DB table name = piatea_<trackName>

Will abort if any of the output files already exist or if track name or table name are already exist
in (grep) trackDb_piatea.ra.
'''

from __future__ import print_function
import argparse
import sys
import re
import os.path
import subprocess

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedFile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='BED infile [default: STDIN]')
    parser.add_argument('-n', '--name', type=str, required=True, 
                        help='name of track, base name of bigBed (<name>.bb), and DB table (piatea_<name>) [required]')
    parser.add_argument('-t', '--trackFolder', type=str, required=True,
                        help='folder from which to publish [required]')
    parser.add_argument('-d', '--desc', type=str, required=False,
                        help='track description [default: name]')
    return parser.parse_args()

    
def exists_file(fileName):
    if not os.path.isfile(fileName): exit("File {} does not exist".format(fileName))

def exists_dir(dirName):
    if not os.path.isdir(dirName): exit("Directory {} does not exist".format(dirName))

def not_exists_file(fileName):
    if os.path.isfile(fileName): exit("File {} already exists".format(fileName))

def process(args):
    exists_dir(args.trackFolder)
    
    sortFileName = "{}.sorted".format(args.bedFile.name)
    bbFileName = "{}.bb".format(args.name)
    for fn in sortFileName, bbFileName:
            not_exists_file(fn)    
    
    with open ('/data/douglas.hoen/trackfile/lyrata/trackDb_piatea.ra') as dbFile:
        with open('/dev/null', 'w') as devnull:
            exists_file(dbFile.name)
            tableName = 'piatea_{}'.format(args.name)
            #tableName = 'piatea'
            checkTableCmd = 'grep {0} {1}'.format(tableName, dbFile.name)
            #print(checkTableCmd)
            tableExists = not subprocess.call(checkTableCmd.split(), stdout=devnull)
            if tableExists:
                exit('{0} already exists in {1}! Exiting.'.format(tableName, dbFile.name))
            
            
            #print(subprocess.check_output(cmd))
        
    
        '''
        sort -k1,1 -k2,2 AL_renamed_scaffolds1to8.bed >AL_renamed_scaffolds1to8_sorted.bed
        bedToBigBed AL_renamed_scaffolds1to8_sorted.bed ../../AL_chrom_sizes.tsv AL_renamed_scaffolds1to8.bb
        cd ~/sharebrowser/tracks/alyrata/permanent/te_search_programs/
        cp /data/douglas.hoen/piatea/tools/repeatmodeler/alyrata/AL_renamed_scaffolds1to8.bb repeatmodeler_repclass.bb
        hgBbiDbLink alyrata piatea_repeatmodeler /data/douglas.hoen/sharebrowser/tracks/alyrata/permanent/te_search_programs/denovo/repeatmodeler_repclass.bb
        '''
    


if __name__ == '__main__':
    args = parse_cl()
    process(args)