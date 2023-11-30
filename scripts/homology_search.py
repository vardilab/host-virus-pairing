import os, sys, subprocess, shlex, re
from collections import defaultdict
import subprocess

hits = []
inputdir = sys.argv[1]
for i in os.listdir(inputdir):
	if i.endswith(".fasta"):
	#if i.endswith("fq.gz"):
		print i
		infile = os.path.join(inputdir, i)
		output = re.sub(".fasta", ".lastout", infile)

		cmd = "lastal -m 100 -f BlastTab -P 32 -u 2 -Q 0 -F 15 /projects/Aylward_Lab/refseq/refSeq207_GVMAGs_appended/refseq207/all_refseq.masked.lastdb "+ infile
		cmd2 = shlex.split(cmd)
		print cmd
		subprocess.call(cmd2, stderr=open("error3.txt", "w"), stdout=open(output, "w"))

		outhandle = open(output, "r")
		qhits = defaultdict(list)
		for i in outhandle.readlines():
			if i.startswith("#"):
				pass
			else:
				line = i.rstrip()
				tabs = line.split("\t")
				query = tabs[0]
				qhits[query].append(query)
				if len(qhits[query]) > 50:
					pass
				else:
					hit = tabs[1]
					###if "|" in hit:
					###	pipe = hit.split("|")
					###	hit = pipe[3]
		
					hits.append(hit)
		outhandle.close()

hit_set = set(hits)
protein2lineage = {}
n=0
print "reading protein2taxa ..."
with open("/projects/Aylward_Lab/refseq/refSeq207_GVMAGs_appended/refseq207/protein2lineage_with_GV_tax.txt", "r") as infile:
#with open("/groups/Aylward_Lab/frank/protein2lineage_with_GV_tax.txt", "r") as infile:
	for i in infile:
		n+=1
		if n > 10000000:
			print "reading protein2taxa ..."
			n = 0

		line = i.rstrip()
		tabs = line.split("\t")
		cluster = tabs[0].rstrip()
		if cluster in hit_set:
			lineage = "\t".join(tabs[1:len(tabs)])
			#print lineage
			#lineage = tabs[0].rstrip()
			protein2lineage[cluster] = lineage
		else:
			pass

print "starting lastout parsing..."

for i in os.listdir(inputdir):
	if i.endswith(".lastout"):
		handle = open(os.path.join(inputdir, i), "r")
		output = re.sub(".lastout", ".lastout.parsed", i)
		outhandle = open(os.path.join(inputdir, output), "w")

		query2hit = defaultdict(list)
		hits = []
		for i in handle.readlines():
			if i.startswith("#"):
				pass
			else:
				line = i.rstrip()
				tabs = line.split("\t")
				query = tabs[0]
		
				hit = tabs[1]
				percid = tabs[2]
				evalue = float(tabs[10])
				bit = float(tabs[11])
				#print query, hit

				###if "|" in hit:
				###	pipe = hit.split("|")
				###	hit = pipe[3]
			
				if len(query2hit[query]) > 0:
					pass

				else:
					try:
						lineage = protein2lineage[hit]
					except:
						lineage = "NA\tNA\tNA\tNA\tNA"

					if evalue > 1e-3:
						pass
					else:
						query2hit[query].append(hit)
						if len(query2hit[query]) > 1:
							pass
						else:
							print len(query2hit[query])					
							#outhandle.write(query +"\t"+ hit +"\t"+ str(evalue) +"\t"+ str(bit) +"\t"+ lineage +"\n")
							outhandle.write(query +"\t"+ hit +"\t"+ str(percid) +"\t"+ str(evalue) +"\t"+ str(bit) +"\t"+ lineage +"\n")
