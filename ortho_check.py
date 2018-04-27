#! /usr/bin/env python3

import csv
import argparse
import Bio.Seq
import Bio.SeqIO
import shmlast.app
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("-o",type=str,required=True,help="Absolute path to orthofuse csv file.")
parser.add_argument("-e",type=float,help="Preferred e value.") #future work
parser.add_argument("-p",type=str,default="uniprot_sprot.fasta.gz",help="Preferred protein database.")
parser.add_argument("-t",type=int,default=1,help="Thread count.") #24 on premise
parser.add_argument("-a",action="store_true",help="In addition to producing output file, tack salvaged TransIDs&sewqs onto orthofuse fasta as well.")
parser.add_argument("-f",type=str,nargs=argparse.REMAINDER,help="List of fasta files.")
args = parser.parse_args()

#do NOT reset every time
genenameseqdb={} #aggregate
genenameseqdbo={} #orthomerged
qualdict={} #aggregate
ortqualdict={} #orthomerged
seqdict={} #aggregate
ortseqdict={} #orthomerged

def builddbs(fasta,seqdicts): #build contig -> seq dictionary for each assembly.fasta file! header = key; value = seq
    handle=gzip.open(fasta,"rt")
    for record in Bio.SeqIO.parse(handle,"fasta"):
        seqdicts[record.id]=record.seq
    handle.close()

#reduce redundancy here too like in translatedict!!!!!!!!!!!
def qualcheck(inputf,whichdict): #iterate through CSV and make dict that takes best hit key=transid:tuple(eval, secondary, KBID) then run THIS through translatedict GET RID OF REDUNDANCY 
    csvname=inputf+".x.uniprot_sprot.fasta.gz.crbl.csv" 
    with open(csvname,"r") as myfile:
        for row in myfile:
            if "bitscore" not in row:
                contig=row.split(",")[8] #NEEDS TO BE TransID aka header; is this right?
                info=(row.split(",")[2],row.split(",")[5],row.split(",")[13].split("|")[2].split("_")[0]) #eval,alignment score (is this ok?),KBID aka gene name?
                if contig not in whichdict: #if gene isn't in there at all yet
                    whichdict[contig]=info #add it
                else: #if the gene name IS already in there, compare values
                    if info[0] < whichdict[contig][0]:  
                        whichdict[contig]=info
                    elif info[0] == whichdict[contig]:
                        if info[1] > whichdict[contig][1]:
                            whichdict[contig]=info                                  


def translatedict(qualdicts,dbdict,retrieve): #use contsigdb to build gene name:(header,seq) dictionary; type is u for unfused or f for fused WHAT GOES IN () TO TAKE EITHER QUALDICT, dbdict
    for eachcontig in qualdicts:
        kbid=qualdicts[eachcontig][2] #KBID aka "gene name"
        if eachcontig == "NODE_10131_length_1841_cov_230.574_g6605_i2":
            print("FOUND IT: {0}".format(kbid))
        #print("'{0}'".format(eachcontig)) 
        getinfo=(eachcontig,retrieve[eachcontig]) #header, seq    
        if kbid not in dbdict.keys():
            dbdict[kbid]=getinfo 
    

#RUNNING THE SCRIPT

for eachone in args.f:
    exec_test = shmlast.app.CRBL(eachone, args.p, n_threads=args.t) #shmlasting unfused fasta files
    exec_test.run(profile_fn=False)
    print(eachone," shmlasted!")    

exec_test = shmlast.app.CRBL(args.o, args.p, n_threads=args.t) #shmlasting orthofused file
exec_test.run(profile_fn=False)
print(args.o," shmlasted!")

for all in args.f:
    builddbs(all,seqdict)

builddbs(args.o,ortseqdict)

#QUALCHECK!!!!!!! !!!!!!!!!!!! !!!!!!!!!!!!! !!!!!!!!
for eachcheck in args.f:
    qualcheck(eachcheck,qualdict) #create qualict (the best hits dictionary for all the unfused files)

qualcheck(args.o,ortqualdict) #create ortqualdict (the best hits directory for all the fused files)


translatedict(qualdict,genenameseqdb,seqdict) #translatedict calls builddbs(each) and then translates it; returns most complete dict after last iteration; for unfused fasta files
translatedict(ortqualdict,genenameseqdbo,ortseqdict) #for orthofused files

with open("nicely_recovered.txt","w") as myfile: #write output file of missing elements
    for every in genenameseqdb:
        if every not in genenameseqdbo:
            myfile.write(">{0}\n{1}\n".format(genenameseqdb[every][0],genenameseqdb[every][1]))

if args.a: #append recovered elements to orthofused file
    with gzip.open(args.o,"a") as ortfile:
        for every in genenameseqdb:
            if every not in genenameseqdbo:
                out_string = ">{0}\n{1}\n".format(str(genenameseqdb[every][0]),str(genenameseqdb[every][1]))
                ortfile.write(out_string.encode("utf-8"))
