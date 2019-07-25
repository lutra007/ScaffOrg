#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 18:14:54 2018

@author: lutra
"""
import re
import sqlite3
import subprocess

def rev_comp(dna):
    trantab = "".maketrans("ACTGactg","TGACtgac")
    return dna[::-1].translate(trantab)

def rotate_plasmid(main,joints,fold,query):
    query_file=main+'query.txt'
    qfile=open(query_file,'w')
    qfile.write(query)
    qfile.close()
    incFIA=fold+'IncFIA.fasta'
    output=main+'output.txt'
    cmd='blastn -query %s -subject %s -out %s -outfmt "6 qstart qend qlen sstart send" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(incFIA,query_file,output)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    result=((open(output).read()).strip()).split('\n')
    result=[r.split('\t') for r in result]
#    print (result)
    result=[r for r in result if r[0]=='1'][0]
#    print ('IncFIA match', result)
    sstart,send=int(result[-2]),int(result[-1])
#    print ('Query len',len(query),end='->')
    if sstart<send:
        query=query[sstart-1:]+query[:sstart-1]
    else:
        print ('FLIP!!!')
        query=query[sstart:]+query[:sstart]
        query=rev_comp(query)
        #reverse joints
        joints=[j+126 for j in joints]
        joints=[len(query)-j+1 for j in joints]
        sstart=len(query)-sstart+1
        
    #process joints
    new_joints=[]
    for j in joints:
        j-=(sstart-1)
        if j<0:
            j+=len(query)
        new_joints.append(j)
#    print (len(query))
    new_joints.sort()
    return query,new_joints

def PrepareGBfromDNA (query):
    text='ORIGIN\n'
    count=1
    fast=[query[i:i+60] for i in range(0, len(query), 60)]
    for f in fast:
        reg=f.strip()
        tab3=' '*(9-len(str(count)))+str(count)+' '
        for i in range(int(len(reg)/10)+1):
            piece=reg[10*i:10*(i+1)]
            tab3=tab3+piece+' '
        tab3+='\n'
        text+=tab3
        count+=60
    return text

def find_coord(main,gene,query,elem,missed):
    query_file=main+'query.txt'
    qfile=open(query_file,'w')
    qfile.write(query)
    qfile.close()
    gene_file=main+'gene.txt'
    gfile=open(gene_file,'w')
    gfile.write(gene)
    gfile.close()
    output=main+'output.txt'
    cmd='blastn -query %s -subject %s -out %s -outfmt "6 qstart qend length qlen sstart send pident" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(gene_file,query_file,output)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    result=(((open(output).read()).strip()).split('\n'))
    if result!=[''] and result:
        result=[[float(rr) for rr in r.split('\t')] for r in result]
        result=[r for r in result if r[2]/r[3]>0.85]
        if len(result)>1:
            result=[r for r in result if r[6]>99.9]
#        if 'hypothetical' in elem:
#            print ('find_coord',elem,len(gene),result)
        matches=[]
        for res in result:
            qstart,qend,align,qlen,sstart,send=[int(r) for r in res[:-1]]
            match='full'
            if align<qlen:
                match='part'
        #        print ('Partial match!')
            nap='ok'
            if sstart>send:
                nap='nope'
                sstart,send=send,sstart
            matches.append([sstart,send,nap,match])
        return matches
    else:
        if missed:
            miss=open(missed,'a')
            entry=elem+'\n'
            print (elem,'is missed')
            miss.write(entry)
            miss.close()
    return []

def CheckProcessedScaffolds(plasmid,reference,main,fold,log_file):
    pl_name=[p for p in ((plasmid.split('/'))[-1]).split('_') if re.search('^p',p)][0]
    nodes={}
    check=[]
    status={}
    duplo={}
    query=(((open(plasmid).read()).strip()).split('>'))[1:]
    output=main+'tmp/output.txt'
    log=open(log_file,'a')
    log.write('\n'+'='*42+'\nCheckProcessedScaffolds\n\n')
    for q in query:
        q=q.split('\n')
        nodes[q[0]]=''.join(q[1:])
        
        query_f=main+'tmp/query.txt'
        que=open(query_f,'w')
        que.write(nodes[q[0]])
        que.close()
        
        cmd='blastn -query %s -subject %s -out %s -outfmt "6 sseqid length qlen qstart qend" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(query_f,plasmid,output)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        result=((open(output).read()).strip()).split('\n')
        for res in result:
            res=res.split('\t')
            node2,length,qlen,qstart,qend=res[0],int(res[1]),int(res[2]),int(res[3]),int(res[4])
            if node2 != q[0]:
                if length/qlen > 0.5:
#                    print (q[0],res)
                    if q[0] not in check:
                        check.append(q[0])
                        status[q[0]]='Duplicate'
                        duplo[q[0]]=[qstart,qend]
    
    cmd='blastn -query %s -subject %s -out %s -outfmt "6 sseqid length qstart qend sstart send qlen" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(reference,plasmid,output)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    result=((open(output).read()).strip()).split('\n')    
    nodes_list=[]
    for res in result:
        res=res.split('\t')
        node=res[0]
        if node not in nodes_list:
            nodes_list.append(node)
            
    for n in nodes:
        if n not in nodes_list:
            if n not in check:
                check.append(n)
                status[n]='Not found'
            else:
                status[n]+=' Not found'
    
    database = fold+'MARA_results.sqlite3'
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    table='Annotation_v4'
    task='SELECT elem,class,nt_seq,strand FROM '+table+' WHERE plasmid="'+pl_name+'"'
    
    annotated={}
    for ch in check:
        annotated[ch]=[]
    
    for row in cur.execute(task):
        el1,cl,nt_seq,nap1=[str(r) for r in row]
        nt_seq=nt_seq.upper()
        
        for ch in check:
            matches=find_coord(main,nt_seq,nodes[ch],el1,False)
            for each in matches:
                st,end,check_nap,match=each
                if el1 not in annotated[ch]:
                    annotated[ch].append(el1)
#                        print (ch,el1,st,end,len(nt_seq),match,end=' ')
                    if ch in duplo:
#                            print ('Duplo!',duplo[ch],st,end)
                        if st in range(duplo[ch][0],duplo[ch][1]+1) and end in range(duplo[ch][0],duplo[ch][1]+1):
#                                print ('Pass!')
                            annotated[ch][-1]+='-Pass'
                        else:
#                                print ('New!')
                            annotated[ch][-1]+='-New'
                    else:
#                            print ('Unique!')
                        annotated[ch][-1]+='-Unique'
                        
    for ch in check:
        if not annotated[ch]:
            status[ch]+=' Not annotated'
        else:
            status[ch]+=(' '+' '.join(annotated[ch]))
        
#    print ('--')
#    for s in status:
#        print (s,status[s])
        
#    print ('--')
    answer=main+'improveSpades/'+(plasmid.split('/'))[-1]
    ans=open(answer,'w')
    for n in nodes:
        add=True
        if n in status:
            if re.search('Not annotated',status[n]) or (not re.search('-New',status[n]) and not re.search('-Unique',status[n])):
                add=False
                log.write(n+' is excluded, this why: '+status[n]+'\n')
                print (n,'is excluded, this why:',status[n])
        if add:
            ans.write('>'+n+'\n'+nodes[n]+'\n')
    ans.close()
    log.write('\nImproved plasmid file: '+answer+'\n\n')
    log.close()
    return answer

def OrganizeScaffolds(plasmid,reference,main,fold,log_file):
#    print (plasmid)
    pl_name=[p for p in ((plasmid.split('/'))[-1]).split('_') if re.search('^p',p)][0]
    ref_name=[p for p in ((reference.split('/'))[-1]).split('_') if re.search('^p',p)][0]
    
    log=open(log_file,'a')
    log.write('\n'+'='*42+'\nOrganizeScaffolds\n')
    
    output=main+'tmp/output.txt'
    cmd='blastn -query %s -subject %s -out %s -outfmt "6 sseqid length qstart qend sstart send qlen" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(reference,plasmid,output)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    result=((open(output).read()).strip()).split('\n')
    nodes=[]
    add_first=0
    nodes_list={}
    for res in result:
#        print (res)
        res=res.split('\t')
        node,align,qst,qend,sst,send,qlen=res[0],int(res[1]),int(res[2]),int(res[3]),int(res[4]),int(res[5]),int(res[6])
        if node not in nodes_list:
            nodes_list[node]='' #to check later if all nodes are added to the final presentation
            add_first=1
        size=(re.findall('_length_([0-9]+)',node))[0]
#        print (node,align,qst,qend,sst,send,qlen)
        if int(size)*0.9<=align or add_first:
            add_first=0
#            print (node)
            if sst<send:
                nodes.append([qst,qend,node,'Forw'])
            else:
                nodes.append([qst,qend,node,'Rev'])
    nodes.sort()
#    for n in nodes:
#        print (n)
    write_nodes=((open(plasmid).read()).strip()).split('>')
    ban=[]
    for wn in write_nodes[1:]:
        wn=wn.split('\n')
        node=wn[0].strip()
        entry=''.join(w.strip() for w in wn[1:])
        if node not in nodes_list:
            nodes_list[node]=''
            add=' '.join(o for o in ['Wow! Abandoned node!', pl_name, node])+'\n'
            ban.append(node)
            log.write(add)
        nodes_list[node]=entry
    final_seq=''
    prst,prend=0,0
    scaff_coords=[]
    log.write('\nReady to organize scaffolds\nQst,qend,node,orientation\n')
    if len(nodes)==1 and nodes[0][1]<0:
        rea=nodes[0][2]
        rearr=nodes_list[rea]
        new=''
        for res in result:
            res=res.split('\t')
            node,align,qst,qend,sst,send,qlen=res[0],int(res[1]),int(res[2]),int(res[3]),int(res[4]),int(res[5]),int(res[6])
            if rea==node:
                if sst>send:
                    send,sst=sst,send
                size=(re.findall('_length_([0-9]+)',node))[0]
                if sst==1:
                    new=new+rearr[-nodes[0][1]:send]
                    ch1=rearr[:-nodes[0][1]]
                elif send==int(size):
                    new=rearr[sst-1:]+new
                    ch2=(rearr[:-nodes[0][1]])
        if ch1==ch2:
            log.write('Perfect match for plasmid'+pl_name+'! Take all info from Ref plasmid '+ref_name+'\n')
#        print (len(new))
        nodes_list[rea]=new
    joints=[]
    for n in nodes:
        add=' '.join(str(o) for o in n)+'\n'
        log.write(add)
        cust,cuend=n[0],n[1]
        node=n[2]
        node_seq=nodes_list[node]
        if n[3]=='Rev':
            node_seq=rev_comp(node_seq)
        if final_seq=='':
            final_seq+=node_seq
            scaff_coords.append([node,1,len(final_seq),0,0,n[3]])    
        elif cust<prend:
            overlap=prend-cust+1
            overlap_reg1=node_seq[:overlap]
            overlap_reg2=final_seq[-overlap:]
            log.write(str(len(overlap_reg1))+','+str(len(overlap_reg2))+'\n')
            if overlap_reg1==overlap_reg2:
                log.write('Overlap!\n')
                final_seq+=node_seq[overlap:]
                joints.append(scaff_coords[-1][2]+1)
                scaff_coords.append([node,scaff_coords[-1][2]+1,len(final_seq),overlap,0,n[3]])
#                print (node,joints[-1])
            else:
                final_seq+=('N'*127)
                final_seq+=node_seq
                joints.append(scaff_coords[-1][2]+1)
                scaff_coords.append([node,scaff_coords[-1][2]+1,len(final_seq),'127N',0,n[3]])
#                print (node,joints[-1])
        else:
            final_seq+=('N'*127)
            final_seq+=node_seq
            joints.append(scaff_coords[-1][2]+1)
            scaff_coords.append([node,scaff_coords[-1][2]+1,len(final_seq),'127N',0,n[3]])
#            print (node,joints[-1])
            
        if cust>cuend:
            print ('FIRST COORDINATE PROBLEM!')
            first_coord=nodes[0][0]
            log.write('First coord '+str(first_coord)+'\n')
            if cuend>first_coord or cuend<0:
                log.write('Close ends!\n')
                print ('CLOSED ENDS!')
                overlap=cuend-first_coord+1
                if cuend<0:
                    overlap=-cuend
                overlap_reg1=final_seq[:overlap]
                overlap_reg2=final_seq[-overlap:]
                if overlap_reg1==overlap_reg2:
                    log.write('Overlap! Initial length '+str(len(final_seq))+'\n')
                    final_seq=final_seq[:-overlap]
                    scaff_coords[-1][2]=len(final_seq)
                    scaff_coords[-1][4]=scaff_coords[-1][2]-scaff_coords[-1][1]+1
        prst,prend=cust,cuend
    last_coord=scaff_coords[-1][2]
#    if ban!=[]:
#        print ('Resubmit',plasmid)
    for banned in ban:
        final_seq+=('N'*127)
        final_seq+=nodes_list[banned]
        plus_coord=len(nodes_list[banned])
        scaff_coords.append([banned,last_coord+1,last_coord+plus_coord,'127N',0,'Forw'])
        joints.append(last_coord+1)
        last_coord+=plus_coord
    if ban:
        joints.append(len(final_seq)+1)
        final_seq+=('N'*127)
    else:
        if len(nodes)>1:
            if final_seq[:127]==final_seq[-127:]:
                log.write('Overlap ends! Well circuled!\n')
                final_seq[:-127]
                joints.append(1)
            else:
                log.write('Not overlapping ends! Add 127N\n')
                joints.append(len(final_seq)+1)
                final_seq+=('N'*127)
        else:
            print ('Single scaffold!')
    log.write('\n')
    add='node, st, end (coords in final_seq), beg_crop, end_crop (if node was cut)\n'
    log.write(add)
    for hist in scaff_coords:
        add='>>>'+pl_name+' '+ref_name+' '+' '.join(str(h) for h in hist)+'\n'
        log.write(add)
#        print (hist)
    name=pl_name+'_'+ref_name
    answer=main+'fasta/'+name+'.fasta'
    ans=open(answer,'w')
    entry='>'+name+'_'+str(len(final_seq))+'bp\n'
    ans.write(entry)
    for i in range(int(len(final_seq)/60)):
        ans.write(final_seq[i*60:(i+1)*60]+'\n')
    ans.write(final_seq[(i+1)*60:])
    ans.close()
    log.close()
    return answer,joints

color_code={'def':'255 255 255'} #white
color_code['hypro']='230 230 230'#light gray
color_code['unchar']='230 230 230'#light gray
color_code['dna']='128 128 128' #gray
color_code['metabol']='128 128 128' #gray
color_code['tRNA']='128 128 128' #gray
color_code['transloc']='128 128 128' #gray

color_code['true_arg']='255 0 0' #red
color_code['antibacterials']='255 102 102' #pink
color_code['antibact']='255 102 102' #pink

color_code['toxin']='0 255 0' #light green
color_code['virdb']='0 77 0' #dark green
color_code['virulence']='0 77 0' #dark green
color_code['vir']='0 77 0' #dark green
color_code['metal']='0 255 255' #aqua

color_code['tra']='0 191 255' #light blue

color_code['Inc']='128 0 128' #purple
color_code['replication']='128 0 128' #purple
color_code['partition']='128 0 128' #purple
color_code['cell_partition']='128 0 128' #purple

color_code['ISel']='128 128 128' #gray
color_code['mobile']='128 128 128' #gray
color_code['IS']='255 140 0' #orange

color_code['phage']='255 204 102' #sand

color_code['joints']='0 0 0' #black
#    color_code['']=''

def GenerateGBfeatures (fasta,gb_file,fold,main,missed):
    fast=((open(fasta).read()).strip()).split('\n')
    length=(re.findall('[0-9]+',fast[0]))[-1]
    query=''.join(fast[1:])
    query=query.upper()
    
    miss=open(missed,'a')
    miss.write('\n\n\n'+'='*42+'\nMissedGenes\n')
    miss.close()
    
    name=(((fasta.split('/'))[-1]).split('.'))[0]
    plasmid=[p for p in name.split('_') if re.search('^p',p)][0]
    if 'CP' in plasmid:
        plasmid=plasmid[1:]
    
    gb=open(gb_file,'w')
    first='LOCUS       '+name+'              '+length+' bp    DNA                         \nFEATURES             Location/Qualifiers\n     source          1..'+length+'\n'
    gb.write(first)
    
    sorty=[]
    
    database = fold+'MARA_results.sqlite3'
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    table='Annotation_v4'
    task='SELECT elem,class,nt_seq,strand,RAST FROM '+table+' WHERE plasmid="'+plasmid+'"'

    checkARGs=fold+'MARA_ARGs.txt'
    args_coll=((open(checkARGs).read()).strip()).split('\n')
    args_db={}
    for a in args_coll:
        aa=a.split('\t')
        args_db[aa[0]]=aa[1:]
    is_library=fold+'IS_library.txt'
    is_lib=((open(is_library).read()).strip()).split('\n')
    isel_db={}
    for il in is_lib:
        iil=il.split('\t')
        isel_db[iil[0]]=iil[1]
    
    rev_nap={'F':'R','R':'F'}
    ready=[]
    for row in cur.execute(task):
        el1,cl,nt_seq,nap1,rast1=[str(r) for r in row]
        nt_seq=nt_seq.upper()
        
        part=False
        for nt in ready:
            if nt_seq in nt:
                part=True
        if nt_seq not in ready and rev_comp(nt_seq) not in ready and not part:
            ready.append(nt_seq)
            ready.append(rev_comp(nt_seq))
            matches=find_coord(main,nt_seq,query,el1,missed)
            real_nap=''
            if nt_seq[1:3]=='TG':
                real_nap='F'
            if nt_seq[-3:-1]=='CA':
                if real_nap:
                    print ('BOTH NAPS?!',el1,nt_seq,'\n')
                real_nap='R'
            if not real_nap:
                print ('NO NAPS?!',el1,nap1,'\n')
                real_nap=nap1
                
            for each in matches:
                el,nap=el1,real_nap
                st,end,check_nap,match=each
                if check_nap!='ok':
                    nap=rev_nap[real_nap]
                if match=='part':
                    if '#' not in el:
                        el+='#'
                
                if len(rast1)>1:
                    el+=('>->->'+rast1)    
                st,end=str(st),str(end)
                sorty.append([int(st),int(end),el,cl,st,end,nap])
    
    sorty.sort()
    exist=[]
    for so in sorty:
        if so not in exist and (exist==[] or (so[0] not in range(exist[-1][0],exist[-1][1]+1) or so[1] not in range(exist[-1][0],exist[-1][1]+1))):
            exist.append(so)
            el,cl,st,end,nap=so[2],so[3],str(so[0]),str(so[1]),so[6]
            ent=st+'..'+end
            
            if nap=='R':
                ent='complement('+ent+')'
            if cl in color_code:
                gene1=el
                if el in isel_db:
                    gene1+=('>->->'+isel_db[el])
                gene1+=('>>>--->>>'+cl)
                entry='     CDS             '+ent+'\n                     /locus_tag="'+gene1+'"\n                     /colour='+color_code[cl]+'\n'
            else:
                entry='     CDS             '+ent+'\n                     /colour='+color_code['def']+'\n'
            gb.write(entry)    
    origin=PrepareGBfromDNA (query)
    gb.write(origin)
    gb.write('//')
    gb.close()
    return int(length)

def GenerateGBfeaturesWJoints (joints,fasta,gb_file,fold,main,missed,flag=False):
    fast=((open(fasta).read()).strip()).split('\n')
    length=(re.findall('[0-9]+',fast[0]))[-1]
    query=''.join(fast[1:])
    query=query.upper()
    
    miss=open(missed,'a')
    miss.write('\n\n\n'+'='*42+'\nMissedGenes\n')
    miss.close()
    
    name=(((fasta.split('/'))[-1]).split('.'))[0]
    plasmid=[p for p in name.split('_') if re.search('^p',p)][0]
    if 'CP' in plasmid:
        plasmid=plasmid[1:]
    
    gb=open(gb_file,'w')
    first='LOCUS       '+name+'              '+length+' bp    DNA                         \nFEATURES             Location/Qualifiers\n     source          1..'+length+'\n'
    gb.write(first)
    
    sorty=[]
    
    database = fold+'MARA_results.sqlite3'
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    table='Annotation_v4'
    task='SELECT elem,class,nt_seq,strand,RAST FROM '+table+' WHERE plasmid="'+plasmid+'"'
    box=['Inc','toxin','true_arg','IS','partition','antibacterials','virulence','virdb','vir','metal']

    checkARGs=fold+'MARA_ARGs.txt'
    args_coll=((open(checkARGs).read()).strip()).split('\n')
    args_db={}
    for a in args_coll:
        aa=a.split('\t')
        args_db[aa[0]]=aa[1:]
    is_library=fold+'IS_library.txt'
    is_lib=((open(is_library).read()).strip()).split('\n')
    isel_db={}
    for il in is_lib:
        iil=il.split('\t')
        isel_db[iil[0]]=iil[1]
    
    rev_nap={'F':'R','R':'F'}
    ready=[]
    for row in cur.execute(task):
        el1,cl,nt_seq,nap1,rast1=[str(r) for r in row]
        nt_seq=nt_seq.upper()
        
        part=False
        for nt in ready:
            if nt_seq in nt:
                part=True
        if nt_seq not in ready and rev_comp(nt_seq) not in ready and not part:
            ready.append(nt_seq)
            ready.append(rev_comp(nt_seq))
            matches=find_coord(main,nt_seq,query,el1,missed)
            real_nap=''
            if nt_seq[1:3]=='TG':
                real_nap='F'
            if nt_seq[-3:-1]=='CA':
                if real_nap:
#                    print ('BOTH NAPS?!',el1,nt_seq,'\n')
                    pass
                real_nap='R'
            if not real_nap:
#                print ('NO NAPS?!',el1,nap1,'\n')
                real_nap=nap1
                
            for each in matches:
                el,nap=el1,real_nap
                st,end,check_nap,match=each
                if check_nap!='ok':
                    nap=rev_nap[real_nap]
                if match=='part':
                    if '#' not in el:
                        el+='#'  
                st,end=str(st),str(end)
                
                if cl in box or el in ('psiA','psiB'):
                    sorty.append([int(st),int(end),el,cl,st,end,nap])
                elif cl in color_code:
                    sorty.append([int(st),int(end),' ',cl,st,end,nap])
                else:
                    print ('New class?!', cl)
    
    sorty.sort()
    exist=[]
    for so in sorty:
        if so not in exist and (exist==[] or (so[0] not in range(exist[-1][0],exist[-1][1]+1) or so[1] not in range(exist[-1][0],exist[-1][1]+1))):
            exist.append(so)
            el,cl,st,end,nap=so[2],so[3],str(so[0]),str(so[1]),so[6]
            ent=st+'..'+end
            
            if nap=='R':
                ent='complement('+ent+')'
            if cl in color_code:
                if cl=='IS' and not re.search('IS26',el) and not re.search('ISEcp1',el) and el not in ('int','tnpR','insB1'):
                    if re.search('#',el):
                        el=(el.split('#'))[0]
                    gene1=isel_db[el]
                else:
                    gene1=el
                entry='     CDS             '+ent+'\n                     /locus_tag="'+gene1+'"\n                     /colour='+color_code[cl]+'\n'
            else:
                entry='     CDS             '+ent+'\n                     /colour='+color_code['def']+'\n'
            gb.write(entry)
            
    if flag:
        for j in joints:
            ent='{}..{}'.format(j,j+126)
            entry='     tRNA            '+ent+'\n                     /colour='+color_code['joints']+'\n'
            gb.write(entry)
    
    origin=PrepareGBfromDNA (query)
    gb.write(origin)
    gb.write('//')
    gb.close()
    return int(length)