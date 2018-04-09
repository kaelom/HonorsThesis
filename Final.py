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


genenameseqdb={} #do NOT reset every time
genenameseqdbo={}


def builddbs(fasta): #build contig -> seq dictionary for each assembly.fasta file! header = key; value = seq
    contsigdb={} #unfused;reset with each file
    print(fasta)
    handle=gzip.open(fasta,"rt")
    print(handle)
    for record in Bio.SeqIO.parse(handle,"fasta"):
        contsigdb[record.id]=record.seq
    handle.close()
    return contsigdb


def translatedict(inputf,type): #use contsigdb to build gene name:(header,seq) dictionary; type is u for unfused or f for fused
    global genenameseqdb
    global genenameseqdbo
    retrieve=builddbs(inputf)
    csvname=inputf+".x.uniprot_sprot.fasta.gz.crbl.csv" #check this
    with open(csvname,"r",newline="") as myfile: #operate on new dict w best hit only (1); make function? 
        for row in myfile:
            if "bitscore" not in row:
                conv=str(row)
                if type=="u":
                    genename=conv.split(",")[13].split("|")[2].split("_")[0]
                    getinfo=(conv.split(",")[8],retrieve[conv.split(",")[8]])
                    if genename not in genenameseqdb:
                        boo=False
                        for each in genenameseqdb.values():
                            if each[0]==getinfo[0]:
                                boo=True
                                break
                        if not boo:
                            genenameseqdb[genename]=getinfo #parsing out gene names to set as key, then using TransID to index contsigdb to set the value as the (header,seq)
#CHECK HERE FOR DUPLICATES 
                else:
                    genename=conv.split(",")[13].split("|")[2].split("_")[0]
       	       	    getinfo=(conv.split(",")[8],retrieve[conv.split(",")[8]])
                    if genename not in genenameseqdbo:
                        boo=False
                        for each in genenameseqdbo.values():
                            if each[0]==getinfo[0]:
                                boo=True
                                break
                        if not boo:
                            genenameseqdbo[genename]=getinfo #parsing out gene names to set as key, then using TransID to index contsigdb to set the value as the (header,seq)
    myfile.close()
    
def convsets(type): #convert dicts to sets
    if type == "a":
        agg=set(genenameseqdb)
        return agg
    else:
        ort=set(genenameseqdbo)
        return ort
#deduplicate
for eachone in args.f:
    exec_test = shmlast.app.CRBL(eachone, args.p, n_threads=args.t) #shmlasting unfused fasta files
    exec_test.run(profile_fn=False)
    print(eachone," shmlasted!")    

exec_test = shmlast.app.CRBL(args.o, args.p, n_threads=args.t) #shmlasting orthofused file
exec_test.run(profile_fn=False)
print(args.o," shmlasted!")

for each in args.f: #for each fasta file that was input
    translatedict(each,"u") #translatedict calls builddbs(each) and then translates it; returns most complete dict after last iteration; for unfused fasta files
    
translatedict(args.o,"f") #for orthofused files

with open("recovered.txt","w") as myfile: #write output file of missing elements
    for every in convsets("a"):
        if every not in convsets("m"):
            myfile.write(">{0}\n{1}\n".format(str(genenameseqdb[every][0]),str(genenameseqdb[every][1])))

if args.a: #append recovered elements to orthofused file
    with gzip.open(args.o,"a") as ortfile:
        for every in convsets("a"):
            if every not in convsets("m"):
                out_string = ">{0}\n{1}\n".format(str(genenameseqdb[every][0]),str(genenameseqdb[every][1]))
                ortfile.write(out_string.encode("utf-8"))
