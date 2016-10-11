#!/usr/bin/env python

from __future__ import division
from datetime import datetime
from subprocess import Popen, PIPE
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.stats as st
import copy,os,glob,sys,random,math,warnings,ast,time,pickle,shutil,operator,scipy
warnings.simplefilter('ignore')
matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
matplotlib.rcParams['font.size'] = 24
matplotlib.style.use('ggplot')


def make_directory(folder,save_data):
	output_folder=(folder+'/')
	if not os.path.isdir(output_folder): 
		os.makedirs(output_folder)
		os.chdir(output_folder)
		os.makedirs("temp")
		if save_data:
			os.makedirs("raw_data")
	else:
		os.chdir(output_folder)
		if not os.path.isdir("temp"): 
			os.makedirs("temp")
		else:
			shutil.rmtree("temp")
			os.makedirs("temp")
		if save_data:
			if not os.path.isdir("raw_data"): 
				os.makedirs("raw_data")
			else:
				shutil.rmtree("raw_data")
				os.makedirs("raw_data")

def get_memory(process):
	import resource
	ram = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	ram = ram/1048576
	print_out(process+' ram : '+str(ram)+' GB') 
def do_with_data(bam1,bam2,depth,lib,case,score_base):
	samtool="/share/apps/samtools/bin/samtools" # Please put your samtools path here
	# samtool, err = Popen(["which","samtools"], stdout=PIPE, stderr=PIPE).communicate()
	# samtool=samtool.replace('\n','')
	def split_chromosome(all_chromosome,bam1,bam2,score_base):
		for chromosome in all_chromosome:
			depth[lib]['1']["+"].update({chromosome:{}})
			depth[lib]['1']["-"].update({chromosome:{}})
			depth[lib]['2']["+"].update({chromosome:{}})
			depth[lib]['2']["-"].update({chromosome:{}})
			print_process("calculate depth "+lib+" : chromosome : "+chromosome+"...")
			data, err=Popen([samtool,"view","../"+bam1,chromosome], stdout=PIPE, stderr=PIPE).communicate()
			data=data.split('\n')
			data=data[:-1]
			depth[lib]['1']["+"][chromosome],depth[lib]['1']["-"][chromosome],total_reads=find_depth(
									data,chromosome,case,score_base)
			depth[lib]['total_reads_1']+=total_reads
			data, err=Popen([samtool,"view","../"+bam2,chromosome], stdout=PIPE, stderr=PIPE).communicate()
			data=data.split('\n')
			data=data[:-1]
			depth[lib]['2']["+"][chromosome],depth[lib]['2']["-"][chromosome],total_reads=find_depth(
									data,chromosome,case,score_base)
			depth[lib]['total_reads_2']+=total_reads
			if memory:
				save_varible([depth[lib]['1']["+"][chromosome],depth[lib]['2']["+"][chromosome]],"temp/"+lib+"/+/depth_"+chromosome)
				depth[lib]['1']["+"][chromosome]={}
				depth[lib]['2']["+"][chromosome]={}
				save_varible([depth[lib]['1']["-"][chromosome],depth[lib]['2']["-"][chromosome]],"temp/"+lib+"/-/depth_"+chromosome)
				depth[lib]['1']["-"][chromosome]={}
				depth[lib]['2']["-"][chromosome]={}		
		return depth
	header,err=Popen([samtool,"view","-H","../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
	header2,err=Popen([samtool,"view","-H","../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
	if err != "" :
		sys.exit(err)
	elif header == "":
		sys.exit("the header of "+bam1+" missing")
	elif header2 == "":
		sys.exit("the header of "+bam2+" missing")
	header=header.split('\n')
	all_chromosome=[]
	for line in header:
		if '@HD' in line:
			issort=line.split("\t")[-1].replace("SO:","")
		elif '@SQ' in line:
			line=line.split('\t')
			chromosome = line[1].replace('SN:','')
			all_chromosome.append(chromosome)
	header2=header2.split('\n')
	for line in header2:
		if '@HD' in line:
			issort2=line.split("\t")[-1].replace("SO:","")
			break
	if issort == "coordinate":
		if not os.path.isfile("../"+bam1+".bai"):
			print_process("create index of " + bam1 +" file...")
			not_use=Popen([samtool,"index","../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
	else:
		print_process("sort "+bam1+" file...")
		not_use=Popen([samtool,"sort","-T","/tmp/"+bam1.replace(".bam","_sort.sorted"),"-o","../"+bam1.replace(".bam","_sort.bam"),"../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
		os.remove("../"+bam1)
		os.rename("../"+bam1.replace(".bam","_sort.bam"),"../"+bam1)
		print_process("create index of " + bam1 +" file...")
		not_use=Popen([samtool,"index","../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
	if issort2 == "coordinate":
		if not os.path.isfile("../"+bam2+".bai"):
			print_process("create index of " + bam2 +" file...")
			not_use=Popen([samtool,"index","../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
	else:
		print_process("sort "+bam1+" file...")
		not_use=Popen([samtool,"sort","-T","/tmp/"+bam2.replace(".bam","_sort.sorted"),"-o","../"+bam2.replace(".bam","_sort.bam"),"../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
		os.remove("../"+bam2)
		os.rename("../"+bam2.replace(".bam","_sort.bam"),"../"+bam2)
		print_process("create index of " + bam2 +" file...")
		not_use=Popen([samtool,"index","../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
	depth = split_chromosome(all_chromosome,bam1,bam2,score_base)
	return depth

def find_depth(input_reads,chromosome,case,score_base):
	position_reads={'+':{},'-':{}}
	position_reads['+'].update({chromosome:{}})
	position_reads['-'].update({chromosome:{}})
	num_position=0
	for line in input_reads:
		inline=line.split('\t')
		try:
			start=int(inline[3])
			chromosome=inline[2]
			CIGAR=inline[5]
			seq=inline[9]
			len_seq=len(seq)
			read_quality=int(inline[4])
			base_score=inline[10]
		except Exception as e:
			print line
			sys.exit(e)
		try:
			FLAG=bin(int(inline[1]))
			FLAG=FLAG[2:]
			if len(FLAG) < 12:
				for x in range(0,12-len(FLAG)):
					FLAG='0'+FLAG
		except:
			FLAG='000000000000'
		if FLAG[-5] == '0' and FLAG[-3] != "1": #strand forward and not unmapped read
			last=0
			num_seq=0
			find_M=False
			for x in range(0,len(CIGAR)):
				if not CIGAR[x].isdigit():	
					if case == 'start':	
						if CIGAR[x] == 'M': #alignment match
							#base=seq[num_seq]
							#score=ord(base_score[num_seq])-33
							#if score >= score_base: 		# check score base
							if True:			# not check score base
								try:
									position_reads['+'][chromosome][start] +=1
									num_position+=1
								except:
									position_reads['+'][chromosome].update({start:1})
									num_position+=1
							break	
						else:
							start=start-int(CIGAR[last:x])
							last=x+1
					elif case == "all":	#whole read
						if CIGAR[x] == 'M':#alignment match
							for position in range(start,start+int(CIGAR[last:x])):
								#base=seq[num_seq]
								#score=ord(base_score[num_seq])-33
								#if score >= score_base: 		# check score base
								if True:			# not check score base
									try:
										position_reads['+'][chromosome][position] +=1
										num_position+=1
									except:
										position_reads['+'][chromosome].update({position:1})
										num_position+=1		
								num_seq+=1
							start=start+int(CIGAR[last:x])
							last=x+1							
						elif CIGAR[x] == 'I': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						else:	
							if CIGAR[x] == 'S':
								num_seq+=int(CIGAR[last:x])
							try:
								start=start+int(CIGAR[last:x])
							except:
								sys.exit('line 93\t'+CIGAR)
							last=x+1
					
					elif case == "end":
						if CIGAR[x] == 'M' :
							start=start+int(CIGAR[last:x])
							last=x+1
							find_M=True
						
						elif CIGAR[x] != 'I':
							if find_M:
								if CIGAR[x] == 'S':
									num_seq+=int(CIGAR[last:x])
								try:
									start=start+int(CIGAR[last:x])
								except:
									sys.exit(CIGAR)
								last=x+1
							else:
								
								last=x+1
						else:
							if find_M:
								num_seq+=int(CIGAR[last:x])
								last=x+1
							else:
								last=x+1	
			if case == 'end':
				#base=seq[num_seq]
				#score=ord(base_score[num_seq])-33
				#if score >= score_base: 		# check score base
				if True:
					try:
						position_reads['+'][chromosome][start-1] +=1
						num_position+=1
					except:
						position_reads['+'][chromosome].update({start-1:1})
						num_position+=1			

		elif FLAG[-5] == '1' and FLAG[-3] != "1": 	#strand reverse and not unmapped read
			find_M=False
			last=0
			num_seq=0
			for x in range(0,len(CIGAR)):
				if not CIGAR[x].isdigit():
					if case == 'start': ## Start_reads
						if CIGAR[x] == 'M' :
							start=start+int(CIGAR[last:x])
							last=x+1
							find_M=True
						
						elif CIGAR[x] != 'I':
							if find_M:
								if CIGAR[x] == 'S':
									num_seq+=int(CIGAR[last:x])
								try:
									start=start+int(CIGAR[last:x])
								except:
									sys.exit(CIGAR)
								last=x+1
							else:
								last=x+1

						else:
							if find_M:
								num_seq+=int(CIGAR[last:x])
								last=x+1
							else:
								last=x+1	
					elif case == 'end':	
						if CIGAR[x] == 'M': #alignment match
							#base=seq[num_seq]
							#score=ord(base_score[num_seq])-33
							#if score >= score_base: 		# check score base
							if True:
								try:
									position_reads['-'][chromosome][start] +=1
									num_position+=1
								except:
									position_reads['-'][chromosome].update({start:1})
									num_position+=1
							break	
						else:
							start=start-int(CIGAR[last:x])
							last=x+1
					elif case == 'all': ## All-reads
						if CIGAR[x] == 'M':
							for position in range(start,start+int(CIGAR[last:x])):
								#base=seq[num_seq]
								#score=ord(base_score[num_seq])-33
								#if score >= score_base:
								if True:				
									try:
										position_reads['-'][chromosome][position] +=1
										num_position+=1
									except:
										position_reads['-'][chromosome].update({position:1})
										num_position+=1		
								num_seq+=1
							start=start+int(CIGAR[last:x])
							last=x+1
							find_M=True		
						elif CIGAR[x] != 'I':
							if find_M:
								if CIGAR[x] == 'S':
									num_seq+=int(CIGAR[last:x])
								try:
									start=start+int(CIGAR[last:x])
								except:
									sys.exit(CIGAR)
								last=x+1
							else:
								last=x+1
						else:
							if find_M:
								num_seq+=int(CIGAR[last:x])
								last=x+1
							else:
								last=x+1
			if case == 'start':
				#base=seq[num_seq]
				#score=ord(base_score[num_seq])-33
				#if score >= score_base: 		# check score base
				if True:
					try:
						position_reads['-'][chromosome][start-1] +=1
						num_position+=1
					except:
						position_reads['-'][chromosome].update({start-1:1})
						num_position+=1
	return [position_reads['+'][chromosome],position_reads['-'][chromosome],num_position]

def check_pseudo(input_depth_1,input_depth_2,reads_fillter,filter_type):
	ratio_all_position={}
	both=0
	only_1=0
	only_2=0
	for strand in input_depth_1:
		ratio_all_position.update({strand:{}})
		for chromosome in sorted(input_depth_1[strand]):
			ratio_all_position[strand].update({chromosome:{}})
			for position in sorted(input_depth_1[strand][chromosome]):
				depth_1=input_depth_1[strand][chromosome][position]
				try:
					depth_2=input_depth_2[strand][chromosome][position]
				except Exception as e:
					if depth_1 >= reads_fillter:
						ratio_all_position[strand][chromosome].update({position:True})
						only_1+=1
				else:
					if filter_type == 'total':
						if depth_1+depth_2 >= reads_fillter:
							ratio_all_position[strand][chromosome].update({position:True})
							both+=1
					else:
						if depth_1 >= reads_fillter and depth_2 >= reads_fillter:
							ratio_all_position[strand][chromosome].update({position:True})
							both+=1
							
	for strand in input_depth_2:
		for chromosome in sorted(input_depth_2[strand]):
			x=set(input_depth_2[strand][chromosome].keys())-set(ratio_all_position[strand][chromosome].keys())
			for position in x:
				if input_depth_2[strand][chromosome][position] >= reads_fillter:
					depth_2=input_depth_2[strand][chromosome][position]
					only_2+=1
	return [both,only_1,only_2]
def find_ratio(input_depth_1,num_1,input_depth_2,num_2,pseudo,lib,save_data,reads_count,filter_type):
	# global memory
	# memory=False
	if save_data:
		write_depth1=open('raw_data/'+lib+"_depth_1.txt",'w')
		write_depth1.write('\t'.join(['strand','chromosome','position','depth'])+'\n')
		write_depth2=open('raw_data/'+lib+"_depth_2.txt",'w')
		write_depth2.write('\t'.join(['strand','chromosome','position','depth'])+'\n')
		write_ratio=open('raw_data/'+lib+"_ratio.txt",'w')
		write_ratio.write('\t'.join(['strand','chromosome','position','ratio'])+'\n')
	ratio_all_position={}
	both=0
	only_1=0
	only_2=0
	num_1_1=num_1/(10**(len(str(num_1))-1))
	num_2_1=num_2/(10**(len(str(num_2))-1))
	num_normalize=float(num_2/num_1)
	for strand in sorted(input_depth_1):
		ratio_all_position.update({strand:{}})
		for chromosome in sorted(input_depth_1[strand]):
			ratio_all_position[strand].update({chromosome:{}})
			if memory:
				input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/depth_"+chromosome)
			for position in sorted(input_depth_1[strand][chromosome]):
				depth_1=input_depth_1[strand][chromosome][position]
				if not pseudo:
					try:
						depth_2=input_depth_2[strand][chromosome][position]
					except Exception as e:
						if depth_1 >= reads_count:
							only_1+=1
							if save_data:
								write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
					else:
						if filter_type == 'total':
							if depth_1+depth_2 >= reads_count:
								ratio_at_position=float(depth_1/depth_2)*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_ratio.write('\t'.join([strand,chromosome,str(position),str(ratio_at_position)])+'\n')
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
						else:
							if depth_1 >= reads_count and depth_2 >= reads_count:
								ratio_at_position=float(depth_1/depth_2)*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_ratio.write('\t'.join([strand,chromosome,str(position),str(ratio_at_position)])+'\n')
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
							
				else:
					try:
						depth_2=input_depth_2[strand][chromosome][position]
					except Exception as e:
						if depth_1 >= reads_count:
							ratio_at_position=float(depth_1+1)*num_normalize
							ratio_all_position[strand][chromosome].update({position:ratio_at_position})
							only_1+=1
							if save_data:
								write_ratio.write('\t'.join([strand,chromosome,str(position),str(ratio_at_position)])+'\n')
								write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
					else:
						if filter_type == 'total':
							if depth_1+depth_2 >= reads_count:
								ratio_at_position=float(depth_1/depth_2)*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_ratio.write('\t'.join([strand,chromosome,str(position),str(ratio_at_position)])+'\n')
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
						else:
							if depth_1 >= reads_count and depth_2 >= reads_count:
								ratio_at_position=float(depth_1/depth_2)*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_ratio.write('\t'.join([strand,chromosome,str(position),str(ratio_at_position)])+'\n')
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
			if memory:
				save_varible([input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]],"temp/"+lib+"/"+strand+"/depth_"+chromosome)
				save_varible(ratio_all_position[strand][chromosome],"temp/"+lib+"/"+strand+"/ratio_"+chromosome)
				input_depth_1[strand][chromosome]={}
				input_depth_2[strand][chromosome]={}
				ratio_all_position[strand][chromosome]={}
	for strand in input_depth_2:
		for chromosome in sorted(input_depth_2[strand]):
			if memory:
				input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/depth_"+chromosome)
				ratio_all_position[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			x=set(input_depth_2[strand][chromosome].keys())-set(ratio_all_position[strand][chromosome].keys())
			for position in x:
				depth_2=input_depth_2[strand][chromosome][position]
				try:
					depth_1=input_depth_1[strand][chromosome][position]
				except:
					if depth_2 >= reads_count:
						only_2+=1
						if save_data:
							write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
						if pseudo:
							ratio_at_position=float(1/(depth_2+1))*num_normalize
							ratio_all_position[strand][chromosome].update({position:ratio_at_position})
							if save_data:
								write_ratio.write('\t'.join([strand,chromosome,str(position),str(ratio_at_position)])+'\n')
			if memory:
				save_varible([input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]],"temp/"+lib+"/"+strand+"/depth_"+chromosome)
				save_varible(ratio_all_position[strand][chromosome],"temp/"+lib+"/"+strand+"/ratio_"+chromosome)
				input_depth_1[strand][chromosome]={}
				input_depth_2[strand][chromosome]={}
				ratio_all_position[strand][chromosome]={}
	if save_data:
		write_depth1.close()
		write_depth2.close()
		write_ratio.close()
	return [ratio_all_position,both,only_1,only_2]
							
def boxcox_transform(enrichment,lamda):
	if lamda > 0:
		enrichment_tranform = (enrichment**lamda - 1) / lamda
	else:
		enrichment_tranform = np.log(enrichment)
	return enrichment_tranform
def boxcox_reverse(enrichment,lamda):
	enrichment_reverse=np.exp(np.log(lamda*enrichment+1)/lamda)
	return enrichment_reverse
def boxcox_data_transform(ratio,lib,p_value,pass_cut_off):
	data=[]
	for strand in ratio:
		for chromosome in ratio[strand]:
			if memory:
				ratio[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			data+=ratio[strand][chromosome].values()
			if memory:
				ratio[strand][chromosome]={}
	takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
	sp=plt.subplot(2, 1,1)
	pp=st.probplot(data, dist=st.norm, plot=sp)
	axes = plt.gca()
	y_min,y_max=axes.get_ylim()
	x_min,x_max=axes.get_xlim()
	r_before=pp[1][2]
	plt.text(0,y_max*0.8,r'$R^2$'+' = '+str(round(r_before**2,4)),fontsize=15)
	plt.title('Before boxcox',fontsize=20)
	sp=plt.subplot(2, 1,2)
	boxcox,lamda=st.boxcox(data)
	pp=st.probplot(boxcox, dist=st.norm, plot=sp)
	x_p=pp[0][0].tolist()
	y_p=pp[0][1].tolist()
	r_after=pp[1][2]
	params = st.norm.fit(boxcox)
	arg = params[:-2]
	loc = params[-2]
	scale = params[-1]
	y=st.norm.ppf(1-p_value, loc=loc, scale=scale, *arg)
	if np.isfinite(y):
		i=y_p.index(takeClosest(y,y_p))
	elif np.isinf(y) and y < 0:
		y=min(boxcox)
		i=1
	elif np.isinf(y) and y > 0:
		y=max(boxcox)
		i=len(boxcox)
	n=len(y_p)
	if i==n:
		val= 0.5**(1/n)
	elif i == 1:
		val = 1 - 0.5**(1/n)
	else:
		val=(i - 0.3175)/(n + 0.365)
	quantiles = st.norm.ppf(val)
	plt.axvline(x=quantiles,color='g',label="p_value "+str(p_value))
	axes = plt.gca()
	y_min,y_max=axes.get_ylim()
	x_min,x_max=axes.get_xlim()
	plt.text(0,y_max*0.8,r'$R^2$'+' = '+str(round(r_after**2,4)),fontsize=15)
	plt.text(-2,y_max*0.6,'lamda = '+str(lamda),fontsize=15)
	plt.text(quantiles,y_max*0.9,"critical value = "+str(round(y,5)),fontsize=15)
	plt.legend(loc='upper left', shadow=True,fontsize=15)
	plt.title('After boxcox',fontsize=15)
	plt.savefig(lib+'_transform.png', dpi=200)
	plt.clf()
	if r_after**2 > pass_cut_off:
		return [lamda,params,y,r_after**2,True]
	else:
		return [lamda,params,y,r_after**2,False]

def transform_all_ratio(ratio,lamda,lib,save_data):
	if save_data:
		write_ratio=open('raw_data/'+lib+"_ratio.txt",'w')
		write_ratio.write('\t'.join(['strand','chromosome','position','ratio(Before boxcox)','ratio(After boxcox)'])+'\n')
	have_to={}
	for strand in sorted(ratio):
		have_to.update({strand:{}})
		for chromosome in sorted(ratio[strand]):
			have_to[strand].update({chromosome:{}})
			if memory:
				ratio[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			for position in sorted(ratio[strand][chromosome]):
				enrichment=ratio[strand][chromosome][position]
				try:
					enrichment_tranform=have_to[strand][chromosome][position]
				except:
					enrichment_tranform=boxcox_transform(enrichment,lamda)
					have_to[strand][chromosome].update({position:enrichment_tranform})
				ratio[strand][chromosome][position]=enrichment_tranform
				if save_data:
					write_ratio.write('\t'.join([strand,chromosome,str(position),str(enrichment),str(enrichment_tranform)])+'\n')
			if memory:
				save_varible(ratio[strand][chromosome],"temp/"+lib+"/"+strand+"/ratio_"+chromosome)
				ratio[strand][chromosome]={}
	if save_data:
		write_ratio.close()
	return ratio
def find_empirical(ratio,percentile,lib):
	data=[]
	for strand in ratio:
		for chromosome in ratio[strand]:
			if memory:
				ratio[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			data+=ratio[strand][chromosome].values()
			if memory:
				ratio[strand][chromosome]={}
	density = st.gaussian_kde(data)
	xs = np.linspace(min(data),max(data),100)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	plt.plot(xs,density(xs),label="empirical distribution")
	plt.fill_between(xs,0, density(xs), facecolor='blue', alpha=0.5)
	peak=st.scoreatpercentile(data,percentile)
	axes = plt.gca()
	y_min,y_max=axes.get_ylim()
	x_min,x_max=axes.get_xlim()
	plt.axvline(x=peak,color='r')
	plt.text(peak,y_max*0.6,"percentile = "+str(percentile)+", critical value = "+str(round(peak,5)),fontsize=15)
	# plt.text(x_max*0.2,y_max*0.8,"total position : "+str(len(data))+" , 5% : "+str(len(data)*0.05),fontsize=15)
	plt.title('Empirical distribution_'+lib)
	plt.xlabel('Ratio')
	plt.ylabel('Probability density')
	plt.legend(loc='upper right', shadow=True)
	plt.savefig(lib+'_empirical.png', dpi=300)
	plt.clf()
	data_unique=list(set(data))
	data_unique.sort()
	all_p_value={}	
	for enrichment in data_unique:
		p_value=st.percentileofscore(data,enrichment)/100
		all_p_value.update({enrichment:p_value})
	return peak,all_p_value

def pie_plot(data,text,name_save,lib):
	def make_autopct(values):
		def my_autopct(pct):
			total = sum(values)
			val = int(round(pct*total/100.0))
			return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
		return my_autopct
	labels = text
	sizes = data
	color = ['salmon','deepskyblue','lightgreen']
	# for x in range(0,len(data)):
	# 	color.append(random.choice(colors.cnames.keys()))
	explode=[0]*len(data)
	explode[0]=0.2
	explode = tuple(explode)  # explode 1st slice
	plt.pie(sizes, labels=labels, colors=color,autopct=make_autopct(sizes),explode=explode, startangle=140)
	plt.axis('equal')
	plt.title(lib)
	plt.savefig(name_save+'.png', dpi=300)
	plt.clf()


def input_gene_gff(input_folder):
	gene_format={'+':{},'-':{}}
	gff_file=glob.glob(os.path.join('../', '*.gff'))
	if len(gff_file) == 0: 
		print input_folder
		sys.exit('no gff')
	elif len(gff_file) > 1:
		sys.exit('there are more than one gff file in this '+input_folder)
	f=open(gff_file[0])

	# gff_file="../gene.gff"
	# if not os.isfile(gff_file):
	# 	sys.exit('gene.gff not contain in '+input_folder)
	# f=open(gff_file)

	f=f.readlines()
	have_chromosome=False
	for line in f:
		line=line.replace('\n','')
		if line[0]=='#':
			if '##FASTA' in line: break	
			if '##sequence-region' in line:
				have_chromosome=True
				inline=line.split()
				len_chromosome=int(inline[3])
				chromosome=inline[1]
				gene_format['+'].update({chromosome:[len_chromosome,{}]})
				gene_format['-'].update({chromosome:[len_chromosome,{}]})
		else:
			if not have_chromosome:
				sys.exit("gff file don't have header line ##sequence-region")
			inline=line.split('\t')
			chromosome=inline[0]
			start=int(inline[3])
			end=int(inline[4])
			size=(abs(end-start))
			strand=inline[6]
			info=inline[-1].split(';')
			gene_name=""
			locus_tag=""
			if inline[2] == 'gene':
				for text in info:
					if 'Name=' in text:
						gene_name=text.replace("Name=",'')
					elif 'locus_tag=' in text:
						locus_tag=text.replace("locus_tag=",'')
				if not gene_name:
					gene_name=str(start)
				if not locus_tag:
					locus_tag="n/a"
				
				gene_format[strand][chromosome][1].update({start:[start,end,size,locus_tag,gene_name]})
				# elif strand =='-':
				# 	gene_format['-'][chromosome][1].update({end:[end,start,size,locus_tag,gene_name]})				
	num_gene=0
	for strand in gene_format:
		for chromosome in gene_format[strand]:
			num_gene+=len(gene_format[strand][chromosome][1])
	return [gene_format,num_gene]

def write_all_position(p_value_cut_off,meta,ratio,depth,total_p_value,params_all,peak_cut_off_all,gene_format,name_out,annotate,empirical_table):
	strand_check=['+','-']
	all_gene_position_sort={}
	if gene_format:
		for strand in gene_format:
			all_gene_position_sort.update({strand:{}})
			for chromosome in gene_format[strand]:
				all_gene_position_sort[strand].update({chromosome:gene_format[strand][chromosome][1].keys()})
				all_gene_position_sort[strand][chromosome].sort()
	write_already={}
	write_file=open(name_out+'_result.txt','w')
	write_gff={}
	if gene_format:
		write_file.write('position\tstrand\tchromosome\tdistance\tregion\tname_gene\tstrand_gene\t')
	else:
		write_file.write('position\tstrand\tchromosome\t')
	for lib in sorted(ratio):
		write_file.write(lib+'(ratio)\t'+lib+'(p-value)\t'+lib+'(depth_1)\t'+lib+'(depth_2)\t*\t')
		if meta != 'combine':
			write_gff[lib]={'+':open(lib+'_forward.gff','w'),'-':open(lib+'_reverse.gff','w')}	
			write_gff[lib]['+'].write("#gff-version 3\n")	
			write_gff[lib]['-'].write("#gff-version 3\n")	
	if meta=="combine":
		write_file.write('combine p-value\n')
		write_gff={'+':open('_'.join(sorted(ratio))+'_combine_p-value_forward.gff','w'),'-':open('_'.join(sorted(ratio))+'_combine_p-value__reverse.gff','w')}
		write_gff['+'].write("#gff-version 3\n")	
		write_gff['-'].write("#gff-version 3\n")	
	elif meta=="replicate":
		write_file.write('total libraries (pass cut-off)\tmax p-value (pass cut-off)\n')
	else:
		write_file.write('total libraries (pass cut-off)\tbest p-value\n')	
	all_lib=sorted(ratio)
	if meta=='combine':
		for strand in total_p_value:
			write_already.update({strand:{}})
			for chromosome in sorted(total_p_value[strand]):
				write_already[strand].update({chromosome:{}})
				if memory:
					for lib2 in all_lib:
						depth[lib2]['1'][strand][chromosome],depth[lib2]['2'][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/depth_"+chromosome)
						ratio[lib2][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/ratio_"+chromosome)
				if gene_format:	
					if annotate == "start" and strand == '-':
						all_position=sorted(total_p_value[strand][chromosome], reverse=True)
					elif annotate == "end" and strand == '+':
						all_position=sorted(total_p_value[strand][chromosome], reverse=True)
					else:
						all_position=sorted(total_p_value[strand][chromosome])
				else:
					all_position=sorted(total_p_value[strand][chromosome])
				for position in all_position:
					combine=total_p_value[strand][chromosome][position]
					try:
						check=write_already[strand][chromosome][position]
					except :
						write_file.write('\t'.join([str(position),strand,chromosome])+'\t')
						write_gff[strand].write('\t'.join([chromosome,name_out,'.',str(position),str(position),'.',strand,'.','p-value='+str(combine)])+'\n')
						if gene_format:
							no_gene=False
							if len(all_gene_position_sort[strand][chromosome]) == 0:
								gene_strand=strand_check[strand_check.index(strand)-1]
								if len(all_gene_position_sort[gene_strand][chromosome]) == 0:
									no_gene=True
							else :
								gene_strand=strand
							if no_gene:
								write_file.write('\t'.join(['n/a','n/a','n/a','n/a'])+'\t')
							else:
								if annotate == "start":
									if gene_strand == '+':
										now=0
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position < start_gene:
												region='Upstream Gene'
												dist=abs(start_gene-position)
												fround=True
											elif position <= end_gene or position == start_gene:
												region='Intra Gene'
												dist=abs(start_gene-position)
												fround=True
											elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
												region='Downstream Gene'
												dist=abs(start_gene-position)
												fround=True
											else:
												now+=1	
									else:
										now=-1
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position > end_gene:
												region='Upstream Gene'
												dist=abs(end_gene-position)
												fround=True
											elif position >= start_gene or position == end_gene:
												region='Intra Gene'
												dist=abs(end_gene-position)
												fround=True
											elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
												region='Downstream Gene'
												dist=abs(end_gene-position)
												fround=True
											else:
												now-=1
								else:
									if gene_strand == '-':
										now=0
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position < start_gene:
												region='Downstream Gene'
												dist=abs(start_gene-position)
												fround=True
											elif position <= end_gene or position == start_gene:
												region='Intra Gene'
												dist=abs(start_gene-position)
												fround=True
											elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
												region='Upstream Gene'
												dist=abs(start_gene-position)
												fround=True
											else:
												now+=1	
									else:
										now=-1
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position > end_gene:
												region='Downstream Gene'
												dist=abs(end_gene-position)
												fround=True
											elif position >= start_gene or position == end_gene:
												region='Intra Gene'
												dist=abs(end_gene-position)
												fround=True
											elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
												region='Upstream Gene'
												dist=abs(end_gene-position)
												fround=True
											else:
												now-=1
								write_file.write('\t'.join([str(dist),region,name_gene,gene_strand])+'\t')
						for lib2 in all_lib:
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except :
								try:
									enrich=depth[lib2]['1'][strand][chromosome][position]
									write_file.write('n/a\tn/a\t'+str(enrich)+'\t0\t*\t')
								except:
									try:
										unenrich=depth[lib2]['2'][strand][chromosome][position]
										write_file.write('n/a\tn/a\t0\t'+str(unenrich)+'\t*\t')
									except:
										write_file.write('n/a\tn/a\tn/a\tn/a\t*\t')
							else:
								try:
									params = params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									p_value=1-empirical_table[lib2][enrichment]
									peak_cut_off=peak_cut_off_all[lib2]
								else:
									peak_cut_off=st.norm.ppf(1-p_value_cut_off,loc=loc,scale=scale, *arg)
									p_value=1-st.norm.cdf(enrichment,loc=loc,scale=scale, *arg)
								try:
									enrich=depth[lib2]['1'][strand][chromosome][position]
								except:
									enrich=0
								try:
									unenrich=depth[lib2]['2'][strand][chromosome][position]
								except:
									unenrich=0
								write_file.write('\t'.join([str(enrichment),str(p_value),str(enrich),str(unenrich)])+"\t*\t")
						write_file.write(str(combine)+'\n')	
						write_already[strand][chromosome].update({position:True})
				if memory:
					for lib2 in all_lib:
						depth[lib2]['1'][strand][chromosome]={}
						depth[lib2]['2'][strand][chromosome]={}
						ratio[lib2][strand][chromosome]={}
	else:
		for strand in total_p_value:
			write_already.update({strand:{}})
			for chromosome in sorted(total_p_value[strand]):
				write_already[strand].update({chromosome:{}})
				if memory:
					for lib2 in all_lib:
						depth[lib2]['1'][strand][chromosome],depth[lib2]['2'][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/depth_"+chromosome)
						ratio[lib2][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/ratio_"+chromosome)
				if gene_format:	
					if annotate == "start" and strand == '-':
						all_position=sorted(total_p_value[strand][chromosome], reverse=True)
					elif annotate == "end" and strand == '+':
						all_position=sorted(total_p_value[strand][chromosome], reverse=True)
					else:
						all_position=sorted(total_p_value[strand][chromosome])
				else:
					all_position=sorted(total_p_value[strand][chromosome])
				for position in all_position:
					try:
						check=write_already[strand][chromosome][position]
					except :
						write_file.write('\t'.join([str(position),strand,chromosome])+'\t')
						if gene_format:
							no_gene=False
							if len(all_gene_position_sort[strand][chromosome]) == 0:
								gene_strand=strand_check[strand_check.index(strand)-1]
								if len(all_gene_position_sort[gene_strand][chromosome]) == 0:
									no_gene=True
							else :
								gene_strand=strand
							if no_gene:
								write_file.write('\t'.join(['n/a','n/a','n/a','n/a'])+'\t')
							else:
								if annotate == "start":
									if gene_strand == '+':
										now=0
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position < start_gene:
												region='Upstream Gene'
												dist=abs(start_gene-position)
												fround=True
											elif position <= end_gene or position == start_gene:
												region='Intra Gene'
												dist=abs(start_gene-position)
												fround=True
											elif position > end_gene and now==len(all_gene_position_sort[gene_strand][chromosome])-1:
												region='Downstream Gene'
												dist=abs(start_gene-position)
												fround=True
											else:
												now+=1	
									else:
										now=-1
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position > end_gene:
												region='Upstream Gene'
												dist=abs(end_gene-position)
												fround=True
											elif position >= start_gene or position == end_gene:
												region='Intra Gene'
												dist=abs(end_gene-position)
												fround=True
											elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
												region='Downstream Gene'
												dist=abs(end_gene-position)
												fround=True
											else:
												now-=1
								else:
									if gene_strand == '-':
										now=0
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position < start_gene:
												region='Downstream Gene'
												dist=abs(start_gene-position)
												fround=True
											elif position <= end_gene or position == start_gene:
												region='Intra Gene'
												dist=abs(start_gene-position)
												fround=True
											elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
												region='Upstream Gene'
												dist=abs(start_gene-position)
												fround=True
											else:
												now+=1	
									else:
										now=-1
										fround=False
										while not fround:
											start_gene=all_gene_position_sort[gene_strand][chromosome][now]
											end_gene=gene_format[gene_strand][chromosome][1][start_gene][1]
											name_gene=gene_format[gene_strand][chromosome][1][start_gene][-1]
											if position > end_gene:
												region='Downstream Gene'
												dist=abs(end_gene-position)
												fround=True
											elif position >= start_gene or position == end_gene:
												region='Intra Gene'
												dist=abs(end_gene-position)
												fround=True
											elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
												region='Upstream Gene'
												dist=abs(end_gene-position)
												fround=True
											else:
												now-=1
								write_file.write('\t'.join([str(dist),region,name_gene,gene_strand])+'\t')
						pass_p=0
						for lib2 in all_lib:
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except :
								try:
									enrich=depth[lib2]['1'][strand][chromosome][position]
									write_file.write('n/a\tn/a\t'+str(enrich)+'\t0\t*\t')
								except:
									try:
										unenrich=depth[lib2]['2'][strand][chromosome][position]
										write_file.write('n/a\tn/a\t0\t'+str(unenrich)+'\t*\t')
									except:
										write_file.write('n/a\tn/a\tn/a\tn/a\t*\t')
							else:
								try:
									params = params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									p_value=1-empirical_table[lib2][enrichment]
									peak_cut_off=peak_cut_off_all[lib2]
								else:
									peak_cut_off=st.norm.ppf(1-p_value_cut_off,loc=loc,scale=scale, *arg)
									p_value=1-st.norm.cdf(enrichment,loc=loc,scale=scale, *arg)
								try:
									enrich=depth[lib2]['1'][strand][chromosome][position]
								except:
									enrich=0
								try:
									unenrich=depth[lib2]['2'][strand][chromosome][position]
								except:
									unenrich=0
								if enrichment >= peak_cut_off:
									pass_p+=1
									write_gff[lib2][strand].write('\t'.join([chromosome,name_out,lib2+"_peak",str(position),str(position),str(enrichment),strand,'.','p-value='+str(p_value),'depth_1='+str(enrich),'depth_2='+str(unenrich)])+'\n')
								write_file.write('\t'.join([str(enrichment),str(p_value),str(enrich),str(unenrich)])+"\t*\t")
						write_file.write(str(pass_p)+'\t'+str(total_p_value[strand][chromosome][position])+'\n')
						write_already[strand][chromosome].update({position:True})
				if memory:
					for lib2 in all_lib:
						depth[lib2]['1'][strand][chromosome]={}
						depth[lib2]['2'][strand][chromosome]={}
						ratio[lib2][strand][chromosome]={}
	write_file.close()

def combine_p_value(ratio,params_all,p_value_cut_off,empirical_table):
	total_p_value={'+':{},'-':{}}
	all_lib=sorted(ratio)
	for strand in ratio[all_lib[0]]:
		for chromosome in ratio[all_lib[0]][strand]:
			total_p_value[strand].update({chromosome:{}})
	for lib in sorted(ratio):
		for strand in ratio[lib]:
			for chromosome in sorted(ratio[lib][strand]):
				if memory:
					for lib2 in all_lib:
						ratio[lib2][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/ratio_"+chromosome)
				for position in sorted(ratio[lib][strand][chromosome]):
					try:
						ckeck=total_p_value[strand][chromosome][position]
					except:
						all_p_value=[]
						for lib2 in sorted(ratio):
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except:
								pass
							else:
								try:
									params=params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									p_value=empirical_table[lib2][enrichment]
								else:
									p_value=st.norm.cdf(enrichment,loc=loc,scale=scale, *arg)
								all_p_value.append(p_value)
						not_use,p_value=st.combine_pvalues(all_p_value)
						p_value=1-p_value
						if p_value <= p_value_cut_off:
							total_p_value[strand][chromosome].update({position:p_value})
				if memory:
					for lib2 in all_lib:
						ratio[lib2][strand][chromosome]={}	
	return total_p_value

def replicate_p_value(ratio,params_all,order,p_value_cut_off,empirical_table):
	all_lib=sorted(ratio)
	total_p_value={'+':{},'-':{}}
	for strand in ratio[all_lib[0]]:
		for chromosome in ratio[all_lib[0]][strand]:
			total_p_value[strand].update({chromosome:{}})
	for lib in sorted(ratio):
		for strand in ratio[lib]:
			for chromosome in sorted(ratio[lib][strand]):
				if memory:
					for lib2 in all_lib:
						ratio[lib2][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/ratio_"+chromosome)
				for position in sorted(ratio[lib][strand][chromosome]):
					try:
						ckeck=total_p_value[strand][chromosome][position]
					except:
						all_p_value=[]
						for lib2 in sorted(ratio):
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except:
								pass
							else:
								try:
									params=params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									p_value=empirical_table[lib2][enrichment]
								else:
									p_value=st.norm.cdf(enrichment,loc=loc,scale=scale, *arg)
								all_p_value.append(p_value)
						all_p_value.sort()
						if order <= len(all_p_value):
							p_value=all_p_value[-order]
							p_value=1-p_value
							if p_value <= p_value_cut_off:
								total_p_value[strand][chromosome].update({position:1-all_p_value[-1]})
				if memory:
					for lib2 in all_lib:
						ratio[lib2][strand][chromosome]={}
	return total_p_value
def union_p_value(ratio,params_all,order,p_value_cut_off,empirical_table):
	all_lib=sorted(ratio)
	total_p_value={'+':{},'-':{}}
	for strand in ratio[all_lib[0]]:
		for chromosome in ratio[all_lib[0]][strand]:
			total_p_value[strand].update({chromosome:{}})
	for lib in sorted(ratio):
		for strand in ratio[lib]:
			for chromosome in sorted(ratio[lib][strand]):
				if memory:
					for lib2 in all_lib:
						ratio[lib2][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/ratio_"+chromosome)
				for position in sorted(ratio[lib][strand][chromosome]):
					try:
						ckeck=total_p_value[strand][chromosome][position]
					except:
						all_p_value=[]
						for lib2 in sorted(ratio):
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except:
								pass
							else:
								try:
									params=params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									p_value=empirical_table[lib2][enrichment]
								else:
									p_value=st.norm.cdf(enrichment,loc=loc,scale=scale, *arg)
								all_p_value.append(p_value)
						all_p_value.sort()
						p_value=all_p_value[-1]
						p_value=1-p_value
						if p_value <= p_value_cut_off:
							total_p_value[strand][chromosome].update({position:1-all_p_value[-1]})
				if memory:
					for lib2 in all_lib:
						ratio[lib2][strand][chromosome]={}
	return total_p_value

def read_per_position(data1,data2,lib):
	density = st.gaussian_kde(data1)
	xs = np.linspace(min(data1),max(data1),100)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	plt.plot(xs,density(xs),label="Library 1")
	plt.fill_between(xs,0, density(xs), facecolor='salmon', alpha=0.5)
	
	density = st.gaussian_kde(data2)
	xs = np.linspace(min(data2),max(data2),100)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	plt.plot(xs,density(xs),label="Library 2")
	plt.fill_between(xs,0, density(xs), facecolor='deepskyblue', alpha=0.5)
	
	plt.title('Reads per Position '+lib)
	plt.xlabel('Reads per Position')
	plt.ylabel('Probability density')
	plt.legend(loc='upper right', shadow=True,fontsize=22)
	plt.savefig(lib+'_Reads_per_Position.png', dpi=300)
	plt.clf()
def print_out(text):
	global last_show
	erase=' '*(len(last_show)+10)
	sys.stdout.write("\033[F")
	print (erase)
	sys.stdout.write("\033[F")
	print (text+"\n\n")
	file_log.write(text+'\n')
	last_show=""
def print_process(text):
	global last_show
	erase=' '*(len(last_show)+10)
	sys.stdout.write("\033[F")
	print (erase)
	sys.stdout.write("\033[F")
	print (text)
	last_show=text
def save_varible(varible,name):
	with open(name+".pickle", 'w') as f:
		pickle.dump(varible, f)

def load_varible(name):
	with open(name+".pickle") as f:
		a = pickle.load(f)
	return a

def main():
	parser = OptionParser()
	parser.add_option("-i", "--input", dest="input", default="", type="string", 
		help="destination directory of input file ***requirement")
	parser.add_option("-r", "--reads", dest="reads", default="all", type="string",
		help="position of reads for calculate (start,all,end) ")
	parser.add_option("-f", "--filter", dest="filter", default="total", type="string",
		help="type of filter reads (total,each) ")
	parser.add_option("-z", "--pseudo", dest="pseudo", default=False, action="store_true", 
		help="pseudo count (default=False)" )
	parser.add_option("-p", "--p_value", dest="p_value", default=0.05, type="float",
		help="p_value cut off (default=0.05)")	
	parser.add_option("-g", "--gene",   dest="gene",default=False, type="string",
		help="add gene annotate to result [start,end] (default=No)" )
	parser.add_option("-m", "--meta", dest="meta", default="union", type="string",
		help="meta analysis of data (default=union) (union,combine,replicate) ")
	parser.add_option("-a", "--replicate", dest="replicate", default=100000, type="int",
		help="number of library pass p-value")
	parser.add_option("-o", "--output", dest="output", default="xxxx", type="string",
		help="output folder name")
	parser.add_option("-q", "--qqplot", dest="pass_normal", default=0.9, type="float",
		help="cut-off r-squared to pass normal distribution after boxcox (default=0.9)")
	parser.add_option("-d", "--distribution", dest="distribution", default="nt", type="string",
		help="type of distribution n(normal distribution) or t(top_rank)")
	parser.add_option("--less_memory", dest="memory", default=False, action="store_true",
		help="use less memory")
	parser.add_option("-c","--reads_count", dest="reads_count", default=1, type="int",
		help="filter reads count")
	parser.add_option("-s","--score_base", dest="score_base", default=1, type="int",
		help="cut-off score-base default=1")
	parser.add_option("--raw_data", dest="raw_data", default=True, action="store_false",
		help="save depth and ratio to text file")

	(options, args) = parser.parse_args()
	
	
	global file_log,memory,last_show,file_time

	
	depth={}
	ratio={}
	params_all={}
	boxcox_lamda={}
	empirical_table={}
	all_lib=[]
	all_file=[]
	percentile_lib={}
	peak_cut_off_all={}
	if not options.input:
		parser.print_help()
		sys.exit()
	input_folder=options.input
	pseudo=options.pseudo
	position_depth=options.reads
	p_value_cut_off=options.p_value
	annotate=options.gene
	meta=options.meta
	replicate=options.replicate
	output=options.output
	pass_normal=options.pass_normal
	distribution=options.distribution
	memory=options.memory
	save_data=options.raw_data
	reads_count=options.reads_count
	score_base=options.score_base
	filter_type=options.filter
	if filter_type == 'each' and pseudo:
		sys.exit('!error : filter type "each" can not use with pseudo count')
	if distribution not in ['nt','n','t']:
		sys.exit('!error : -d must be "n" or "t"')
	if meta not in ['combine','replicate','union']:
		sys.exit('!error : -m must be "combine", "replicate" or "union"')
	if annotate and annotate not in ['start', 'end']:
		sys.exit('!error : -g must be "start" or "end"')
	if position_depth not in ['start','all','end']:
		sys.exit('!error : -r must be "start", "all" or "end"')
	else:
		if input_folder[-1] != '/':
			input_folder=input_folder+'/'
		if os.path.isdir(input_folder):
			all_file+=glob.glob(os.path.join(input_folder, '*.bam'))
		else:
			sys.exit('!error : No such directory ('+input_folder+')')
	all_file.sort()
	for x in range(len(all_file)):
		all_file[x]=all_file[x].split("/")[-1]

	if len(all_file) == 0:
		sys.exit("!error : No such sam/bam file in this folder("+input_folder+")")
	else:
		file1=[]
		file2=[]
		for file in all_file:
			if '_1.' in file or '_2.' in file:
				if '_1.' in file:
					file1.append(file[:-6])
				elif '_2.' in file:
					file2.append(file[:-6])
			else :
				sys.exit('!error : Invalid file name ('+file+'), File name must end with "_1 or _2" (exp. ecoli_1.bam or ecoli_2.bam)')
		if len(set(file1)-set(file2)) != 0:
			for file in set(file1)-set(file2):
				print("No such file "+file+"_2")
			sys.exit()
		elif len(set(file2)-set(file1)) != 0:
			for file in set(file2)-set(file1):
				print("No such file "+file+"_1")
			sys.exit()
	
	if meta == 'replicate' and len(all_file)/2 < replicate and replicate != 100000:
		sys.exit('!error : number of library pass p-value more than library')
	elif meta == 'replicate' and replicate < 1:
		sys.exit('!error : number of library pass p-value must more than 0')
	if meta == 'replicate' and replicate == 100000:
		replicate=int(len(all_file)/2)
	if pass_normal > 1:
		pass_normal=1
	if len(all_file) == 2:
		meta='union'
	# ********************start here****************
	time_start=str(datetime.now())[:-7].replace(' ','_').replace(':','_')
	
	if output != 'xxxx':
		make_directory(input_folder+output,save_data)
	else:
		make_directory(input_folder+output+'_'+time_start,save_data)
	time_start=time.time()	
	file_log=open('log.txt','w')
	file_time=open("time.txt",'w')
	last_show=''	
	print_out("\n\n--------------------------")
	print_out("|"+(str(datetime.now())[:-7])+"|")
	print_out("\n\n--------------------------")
	print_out("input_folder : "+input_folder)
	print_out("output folder : " +os.getcwd())
	print_out("There are " +str(int(len(all_file)/2))+ " libraries")
	for file in file1:
		print_out('\t'+file)
	if len(all_file)/2 > 1:
		print_out("meta analysis : "+meta)
	print_out("p-value cut off : "+ str(p_value_cut_off))
	print_out("position of reads : "+position_depth)
	if distribution == "nt":
		print_out("distribution : normal distribution and top rank")
	elif distribution == 'n':
		print_out("distribution : normal distribution")
	else:
		print_out("distribution : top rank")
	print_out("r-squared to pass normal distribution after boxcox : "+str(pass_normal))
	if meta == "replicate":
		print_out("total libraries pass p-value : "+str(replicate))
	print_out("type of filter : "+str(filter_type))
	print_out("number of reads count : "+str(reads_count))
	if pseudo: 
		print_out("pseudo count : Yes")
	else: 
		print_out("pseudo count : No")
	if annotate : 
		print_out("put annotate gene : "+annotate)
	else: 
		print_out("put annotate gene : No")
	if memory:
		print_out("less memory : Yes (use long time)")
	else:
		print_out("less memory : NO")
	print_out("\n****************************\n\n")
	if annotate:
		print_process("import gene annotate...")
		gene_format,num_gene=input_gene_gff(input_folder)

	for file1 in all_file:
		if '_1.' in file1.split('/')[-1]:
			lib=file1[:-6]
			if lib not in all_lib: 
				all_lib.append(lib)
				depth[lib]={'1':{"+":{},"-":{}},'2':{"+":{},"-":{}},'total_reads_1':0,'total_reads_2':0}
				ratio[lib]={}
				file2=file1.replace('_1.','_2.')
				if memory:
					if not os.path.isdir("temp/"+lib):
						os.makedirs("temp/"+lib)
					if not os.path.isdir("temp/"+lib+"/+"):
						os.makedirs("temp/"+lib+"/+")
						os.makedirs("temp/"+lib+"/-")
				depth = do_with_data(file1,file2,depth,lib,position_depth,score_base)
				print_process("calculate ratio "+lib+"...")
				if pseudo:
					both,only_1,only_2=check_pseudo(depth[lib]['1'],depth[lib]['2'],reads_count,filter_type)
					depth[lib]['total_reads_1']+=only_2
					depth[lib]['total_reads_2']+=only_1
				
				ratio[lib],both,only_1,only_2=find_ratio(
									depth[lib]['1'],depth[lib]['total_reads_1'],
									depth[lib]['2'],depth[lib]['total_reads_2'],
									pseudo,
									lib,
									save_data,
									reads_count,
									filter_type
									)
				pie_plot([both,only_1,only_2],['both','only_1','only_2'],lib+'_type_of_reads',lib)
				if distribution == "nt":
					print_process("transform data with boxcox method "+ lib +"...")
					boxcox_lamda[lib],params_all[lib],peak_cut_off_all[lib],r2,normal = boxcox_data_transform(ratio[lib],lib,p_value_cut_off,pass_normal)
					print_out(lib+" : boxcox r-squared = "+str(r2))
					if normal: 
						mu, std = params_all[lib]
						print_out(lib+" : mu = "+str(mu))
						print_out(lib+" : std = "+str(std))
						print_out(lib+" : lamda = "+str(boxcox_lamda[lib]))
						print_out(lib+" : peak cut-off = "+str(peak_cut_off_all[lib]))
						ratio[lib] = transform_all_ratio(ratio[lib],boxcox_lamda[lib],lib,save_data)
					else : 
						del params_all[lib],boxcox_lamda[lib]
						peak_cut_off_all[lib],empirical_table[lib]=find_empirical(ratio[lib],(1-p_value_cut_off)*100,lib)
						print_out(lib+" : peak cut-off = "+str(peak_cut_off_all[lib]))
				elif distribution == "n":
					print_process("transform data with boxcox method "+ lib +"...")
					pass_normal=0
					boxcox_lamda[lib],params_all[lib],peak_cut_off_all[lib],r2,normal = boxcox_data_transform(ratio[lib],lib,p_value_cut_off,pass_normal)
					mu, std = params_all[lib]
					print_out(lib+" : boxcox r-squared = "+str(r2))
					print_out(lib+" : mu = "+str(mu))
					print_out(lib+" : std = "+str(std))
					print_out(lib+" : lamda = "+str(boxcox_lamda[lib]))
					print_out(lib+" : peak cut-off = "+str(peak_cut_off_all[lib]))
					ratio[lib] = transform_all_ratio(ratio[lib],boxcox_lamda[lib],lib,save_data)
				else:
					peak_cut_off_all[lib],empirical_table[lib]=find_empirical(ratio[lib],(1-p_value_cut_off)*100,lib)
					print_out(lib+" : peak cut-off = "+str(peak_cut_off_all[lib]))
	total_p_value={}
	save_varible(ratio,"temp/ratio")
	save_varible(depth,"temp/depth")
	save_varible(params_all,"temp/params_all")
	if meta == 'combine' :
		print_process("combine p-value...")
		total_p_value=combine_p_value(ratio,params_all,p_value_cut_off,empirical_table)
		# save_varible(total_p_value,"temp/total_p_value")
	elif meta == 'replicate':
		print_process("replicate libraries..")
		total_p_value=replicate_p_value(ratio,params_all,replicate,p_value_cut_off,empirical_table)
	else:
		print_process("union libraries..")
		total_p_value=union_p_value(ratio,params_all,replicate,p_value_cut_off,empirical_table)
	if annotate:
		print_process("write result...")
		write_all_position(p_value_cut_off,meta,ratio,depth,total_p_value,params_all,peak_cut_off_all,gene_format,output,annotate,empirical_table)
	else:
		print_process("write result...")
		write_all_position(p_value_cut_off,meta,ratio,depth,total_p_value,params_all,peak_cut_off_all,False,output,annotate,empirical_table)
	
	print_process("..Done..                            ")
	print_out("\n--------------------------")
	print_out("|"+(str(datetime.now())[:-7])+"|")
	print_out("\n--------------------------")

	time_use=(time.time()-time_start)/60
	get_memory("all")
	file_time.write("all time of program : "+str(time_use)+'min \n')
	file_time.close()
	file_log.close()
	
if __name__ == '__main__' :
		main()

