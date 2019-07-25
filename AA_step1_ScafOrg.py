#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 17:52:35 2018

@author: lutra
"""
import re
import glob
import subprocess
from AA_step1_functions import *

group='U1_gr'
fold='/home/shiri/plasmid_project/Current_work/ScafOrg_program/'
main=fold+group+'/'

for mkdir in ('improveSpades','tmp','fasta','rotated_fasta','logs','genbank','genbank_wjoints','easyfig'):
    cmd='mkdir {}{}'.format(main,mkdir)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()

spades=main+'spades/'

ref=glob.glob(main+'*')

choose_ref={'gr1':'ref','gr2':'ref2','U1_gr':'ref'}
ref=[r for r in ref if choose_ref[group] in r][0]
ref_name=[p for p in ((ref.split('/'))[-1]).split('_') if re.search('^p',p)][0]
print (ref_name)

res={}
for pl in glob.glob(spades+'*')[:]:
    pl_name=[p for p in ((pl.split('/'))[-1]).split('_') if re.search('^p',p)][0]
    print (pl_name)
    
    file_name='{}_{}'.format(pl_name,ref_name)
    if True:
        log_file=main+'logs/'+file_name+'.log'
        log=open(log_file,'w')
        entry='Organize: {}\nReference: {}\n\nPlasmid file: {}\nReference file: {}\n\n'.format(pl_name,ref_name,pl,ref)
        log.write(entry)
        log.close()
        
        imp_pl=CheckProcessedScaffolds(pl,ref,main,fold,log_file)
    #    print (imp_pl)
    #    imp_pl=pl
        fasta,joints=OrganizeScaffolds(imp_pl,ref,main,fold,log_file)
    #    print (fasta)
        
        query=''.join((((open(fasta).read()).strip()).split('\n'))[1:])
        print (pl_name,len(query))
        
        query,joints=rotate_plasmid(main,joints,fold,query)
        
        que_file='{}rotated_fasta/{}.fasta'.format(main,file_name)
        que=open(que_file,'w')
        label='>{}_length_{}bp\n'.format(file_name,len(query))
        que.write(label)
        for i in range(int(len(query)/60)):
            que.write(query[i*60:(i+1)*60]+'\n')
        que.write(query[(i+1)*60:])
        que.close()
        
    #    for j in joints:
    #        print (j,query[j-1:j+126],len(query[j-1:j+126]))
        
        gb_file='{}genbank/{}.gb'.format(main,file_name)
#        GenerateGBfeatures(que_file,gb_file,fold,main,log_file)
        gb_file='{}genbank_wjoints/{}.gb'.format(main,file_name)
        
        #last flag for scaffolds separation: TRUE - show borders, False - hide them
        GenerateGBfeaturesWJoints(joints,que_file,gb_file,fold,main,log_file,False)
    gb_file='{}genbank_wjoints/{}.gb'.format(main,file_name)
    res[pl_name]=gb_file

orders={'gr1':['pCP013657','p461','pB316','p396','p418','pU60','p76','p113','pB28'],'gr2':['pCP021180','p425','p131','pU10','pU23','pB9','p26','pB14','pU17','pU1','p56'],'U1_gr':['pU2']}
cmd='python2.7 /home/shiri/plasmid_project/tools/Easyfig-lutra/Easyfig.py -o {}all_{}.svg -svg -width 8821 -ann_height 250 -blast_height 900 -f1 T -f2 10000 -min_length 1000 -aln left -blast_col 0 191 255 0 0 255 -blast_col_inv 255 215 0 255 140 0 -bo F -f CDS 0 0 0 rect -f tRNA 0 0 0 pointer -glt 5 -genet 20 -legend both -leg_name locus_tag -uncomp T'.format(main,ref_name)
#for p in ['pCP013657','p461','pB316','p396','p418','pU60','p76','p113','pB28']:
figrefs={'gr1':'pCP013657','gr2':'p131','U1_gr':'pU1'}
for p in orders[group]:
    add=res[p]
    cmd+=(' '+add)
    figref=figrefs[group]
    gb_ref=res[figref]
    cmd1='python2.7 /home/shiri/plasmid_project/tools/Easyfig-lutra/Easyfig.py -o {}easyfig/{}_{}.svg -svg -width 8821 -ann_height 250 -blast_height 900 -f1 T -f2 10000 -min_length 1000 -aln left -blast_col 0 191 255 0 0 255 -blast_col_inv 255 215 0 255 140 0 -bo F -f CDS 0 0 0 rect -f tRNA 0 0 0 pointer -glt 5 -genet 20 -legend both -leg_name locus_tag -uncomp T {} {}'.format(main,p,figref,add,gb_ref)
    process = subprocess.Popen(cmd1, shell=True)
    process.wait()
    
#print (cmd)
process = subprocess.Popen(cmd, shell=True)
process.wait()
