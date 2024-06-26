#!/usr/bin/env/python
import sys, os, math, optparse, glob, operator, random
from optparse import Values
import pandas as pd

options = None
model_df = None

#########################################################################################################################

class rec_entry:
    def __init__(self,chromosome,begin_position,end_position,combined_rate,genetic_map_dist):
        self.chromosome=chromosome
        self.begin_position=int(begin_position)
        self.end_position=int(end_position)
        self.combined_rate=float(combined_rate)
        self.genetic_map_dist=float(genetic_map_dist)

class model_class:
    def __init__(self,ancestors,ml,meioses,ml_segment_divider,total_segment_count,segments,pair_id):
        self.ancestors=int(ancestors)
        self.ml=float(ml)
        self.meioses=int(meioses)
        self.ml_segment_divider=int(ml_segment_divider)
        self.total_segment_count=int(total_segment_count)
        self.segments=segments
        self.pair_id=pair_id

def get_chromosome(filename):
    chromosome_tmp=filename.split("chr")
    if chromosome_tmp[1][1] in ('0','1','2','3','4','5','6','7','8','9'):
        chromosome='chr'+chromosome_tmp[1][0:2]
    elif chromosome_tmp[1][1].isalpha():
        chromosome='chr'+chromosome_tmp[1][0:2]
    else:
        chromosome='chr'+chromosome_tmp[1][0]
    return chromosome

def factln(x): 
    if x==0:
        return 0
    else:
        return math.log(x) + factln(x-1)

def get_overlap(first_segment,second_segment):
    if first_segment[0]<=second_segment[1] and first_segment[1]>=second_segment[0]:
        if first_segment[0]<=second_segment[0]:
            overlap_begin=second_segment[0]
        else:
            overlap_begin=first_segment[0]
        if first_segment[1]>=second_segment[1]:
            overlap_end=second_segment[1]
        else:
            overlap_end=first_segment[1]
        v_return=first_segment[2]*(overlap_end-overlap_begin)/(first_segment[1]-first_segment[0])
    else:
        v_return=0.0
    return v_return

def get_total_overlap(i,segments,expected_length=-1):
    global options
    first_segment=segments[i]
    observed_length=0.0
    j=0
    while j<len(segments) and segments[j][0]<=first_segment[1]:
        if i!=j:
            observed_length+=get_overlap(segments[i],segments[j])
        if expected_length>0 and observed_length/expected_length>options.mask_region_threshold:
            break
        j+=1
    return observed_length

def simulate_segments(segments):
    global options
    iterations=options.max_region_simulation_count
    for iter in range(iterations):
        if options.verbose:
            print ("...iteration "+str(iter)+" of "+str(iterations))
        for i in range(len(segments)):
            segment=segments[i]
            segment[0]=float(random.random()*(options.rec_per_meioses*100-segment[2]))
            segment[1]=float(segment[0]+segment[2])
        segments.sort(key=operator.itemgetter(0))
        for i in range(len(segments)):
            segment=segments[i]
            observed_length=get_total_overlap(i,segments)
            if observed_length>=segment[3]:
                segment[4]+=1
            segment[8]+=observed_length/float(iterations)

def gammln(xx):
    tmp=ser=0.0
    cof=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179E-2,-0.5395239384953E-5]
    y=x=float(xx)
    tmp=x+5.5
    tmp-=(x+0.5)*math.log(tmp)
    ser=1.000000000190015
    for j in range(6):
        y+=1
        ser+=cof[j]/y
    return -tmp+math.log(2.5066282746310005*ser/x)

def gser(a,x):
    ITMAX=101
    EPS=3.0E-7
    v_return=[]
    gln=gammln(float(a))
    if x<0:
        raise RuntimeError("%prog: x less than 0")
    elif x==0:
        return [0.0]
    else:
        ap=float(a)
        v_del=sum=1.0/a
        for n in range(1,ITMAX):
            ap+=1
            v_del*=x/ap            
            sum+=v_del
            if abs(v_del)<abs(sum)*EPS:
                v_return=[sum*math.exp(-x+a*math.log(x)-gln),gln] 
                break
        if v_return==[]:
            raise RuntimeError("%prog: a too large, ITMAX too small")
        else:
            return v_return

def gcf(a,x):
    ITMAX=101
    EPS=3.0E-7
    FPMIN=1.0E-30
    gln=gammln(a)
    b=x+1.0-a
    c=1.0/FPMIN
    d=1.0/b
    h=d
    for i in range(1,ITMAX):
        an=-i*(i-a)
        b+=2.0
        d=an*d+b
        if abs(d)<FPMIN:
            d=FPMIN
        c=b+an/c
        if abs(c)<FPMIN:
            c=FPMIN
        d=1.0/d
        v_del=d*c
        h*=v_del
        if abs(v_del-1.0)<EPS:
            break
    if i==ITMAX-1:
        raise RuntimeError("%prog: a too large, ITMAX too small")
    gammcf=math.exp(-x+a*math.log(x)-gln)*h
    return [gammcf,gln] 
    
def gammp(a,x):
    if x>=0 and a>0:
        if x<a+1:
            return gser(a,x)[0]
        else:
            return 1.0-gcf(a,x)[0]

def chiinv(a,df,digits=10):
    low=1*df
    p=gammp(df/2.0,low)
    while p>a:
        low=low*0.5
        p=gammp(df/2.0,low)
    high=2*df
    p=gammp(df/2.0,high)
    while p<a:
        high=high*2
        p=gammp(df/2.0,high)
    while round(low,digits)!=round(high,digits):
        current=(low+high)/2.0
        p=gammp(df/2.0,current)
        if p<a:
            low=current
        else:
            high=current
    return 2*current

def get_emp_shared_segment_ll(emp_segment_lambda,segment_count):
    if emp_segment_lambda>0:
        return -emp_segment_lambda+segment_count*math.log(emp_segment_lambda)-factln(segment_count)
    else:
        return min_ll_constant

def set_confidence(conf_level):
    cstat=chiinv(conf_level,2)/2.0
    return [cstat,conf_level]

def background_ll(segment_list,emp_segment_lambda,emp_lambda,ascertained_segments=[],related_segments=[],related_asc_segments=[]):
    global options
    bn=0
    for segment in segment_list:
        if segment>0:
            bn+=1
    if options.adjust_pop_dist=="false":
        b_ll=get_emp_shared_segment_ll(emp_segment_lambda,bn)
    else:
        related_sum=0.0
        for segment in related_segments:
            related_sum+=segment
        for segment in related_asc_segments:
            related_sum+=segment
        unrelated_prop=1-related_sum/(genetic_map)
        b_ll=get_emp_shared_segment_ll(unrelated_prop*emp_segment_lambda,bn)
    for segment in segment_list:
        segment_size=abs(segment)-options.min_cm 
        b_ll+=math.log(emp_lambda)-segment_size*emp_lambda
    for segment in ascertained_segments:
        segment_size=abs(segment)-options.min_cm 
        b_ll+=math.log(emp_lambda)-segment_size*emp_lambda
    return b_ll

def related_0p_ll(segment_list,n,ascertained_segments=[]):
    global options
    rn=0
    for segment in segment_list:
        rn+=1
    bin_n=(genetic_map/100.0)*(n-1)+options.number_of_chromosomes
    if n>1:
        #expected_length=(100.0/(n-1))
        expected_length=genetic_map/bin_n
    else:
        expected_length=100*options.rec_per_meioses/float(options.number_of_chromosomes)
    rc=math.log(expected_length)
    p=math.exp(-options.min_cm/(expected_length))*0.5**(n-1)
    lam=bin_n*p
    r_ll=-lam+rn*math.log(lam)-factln(rn)
    for segment in segment_list:
        segment_size=abs(segment)-options.min_cm
        r_ll+=-rc-segment_size/expected_length
    for segment in ascertained_segments:
        segment_size=abs(segment)-options.min_cm
        r_ll+=-rc-segment_size/expected_length
    return r_ll

def related_1p_ll(segment_list,n,ascertained_segments=[],genome_proportion=1.0):
    global options
    rn=0
    for segment in segment_list:
        rn+=1
    bin_n=(genetic_map/100.0)*n+options.number_of_chromosomes
    #expected_length=(100.0/n)
    expected_length=genetic_map/bin_n
    rc=math.log(expected_length)
    p=math.exp(-options.min_cm/(expected_length))*0.5**(n-1)
    lam=bin_n*p*genome_proportion
    r_ll=-lam+rn*math.log(lam)-factln(rn)
    for segment in segment_list:
        segment_size=abs(segment)-options.min_cm
        r_ll+=-rc-segment_size/expected_length
    for segment in ascertained_segments:
        segment_size=abs(segment)-options.min_cm
        r_ll+=-rc-segment_size/expected_length
    if n>1:
        return r_ll
    else:
        return min_ll_constant

def related_2p_ll(segment_list,n,ascertained_segments=[],genome_proportion=1.0):
    global options
    rn=0
    for segment in segment_list:
        rn+=1
    if n in (2,3):
        p_inh=0.75
        lam=p_inh*options.number_of_chromosomes+((genetic_map/100.0)*2*2)*p_inh*(1-p_inh)
        if n==2:
            expected_map=genetic_map*0.75
        else:
            expected_map=genetic_map*0.5
        expected_length=expected_map/lam
    else:
        bin_n=2*((genetic_map/100.0)*n+options.number_of_chromosomes)
        expected_length=100.0/n #2*genetic_map/bin_n
        p=math.exp(-options.min_cm*(1/expected_length))*0.5**(n-1)
        lam=bin_n*p*genome_proportion
    rc=math.log(expected_length)
    l=2.0/(expected_length)
    r_ll=-lam+rn*math.log(lam)-factln(rn)
    for segment in segment_list:
            segment_size=abs(segment)-options.min_cm
            if n==2:
                l_sum=0.0
                k_fact=1.0
                for k in range(1,10):
                    k_fact=k_fact*float(max(k-1,1))
                    l_sum+=0.5**k*segment_size**(k-1)*math.exp(-segment_size*l)*l**k*(1/k_fact)
                r_ll+=math.log(l_sum)
            else:
                r_ll+=-rc-segment_size/expected_length
    for segment in ascertained_segments:
            segment_size=abs(segment)-options.min_cm
            if n==2:
                l_sum=0.0
                k_fact=1.0
                for k in range(1,10):
                    k_fact=k_fact*float(max(k-1,1))
                    l_sum+=0.5**k*segment_size**(k-1)*math.exp(-segment_size*l)*l**k*(1/k_fact)
                r_ll+=math.log(l_sum)
            else:
                r_ll+=-rc-segment_size/expected_length
    if n==1:
        return min_ll_constant
    else:
        return r_ll

def ibd2_sib_ll(segment_list):
    global options
    rn=0
    k_count=0
    for segment in segment_list:
        rn+=1
    rc=math.log(25.0)
    bin_n=((genetic_map/100.0)+options.number_of_chromosomes)/4.0
    p=math.exp(-options.min_cm/25.0)
    lam=bin_n*p
    r_ll=-lam+rn*math.log(lam)-factln(rn)
    for segment in segment_list:
        segment_size=abs(segment)-options.min_cm
        r_ll+=-rc-segment_size/(25.0)
    return r_ll

def add_segment(ind_sharing,ind_id,cm,controls="no"):
    global options
    if abs(cm)>=options.min_cm:
        if controls=="no" or cm<=options.max_cm:
            if ind_id in ind_sharing: 
                ind_sharing[ind_id].append(cm)
            else:
                ind_sharing[ind_id]=[cm]

def get_masked_coordinates(chromosome,begin_position,end_position,masked_segments_dict):
    global options
    new_begin=begin_position
    new_end=end_position
    if chromosome in masked_segments_dict:
        for segment in masked_segments_dict[chromosome]:
            if new_begin<=segment[1] and new_end>=segment[0]:
                if new_begin>segment[0]-options.mask_region_cross_length:
                    new_begin=min(segment[1],new_end)
                if new_end<segment[1]+options.mask_region_cross_length:
                    new_end=max(segment[0],new_begin)
    return [new_begin,new_end]

def process_segment(chromosome,ascertained_dict,sharing_dict,ibd2_dict,ind_id,cm,controls,begin_position,end_position,recombination_rates,IBD2,control_segments,masked_segments_dict,masked_sum):
    global options
    if ind_id not in masked_sum:
        masked_sum[ind_id]=0.0
    if chromosome==options.ascertained_chromosome and options.ascertained_position>=begin_position and options.ascertained_position<=end_position and controls=="no" and IBD2=="no":
        [new_begin,new_end]=get_masked_coordinates(chromosome,begin_position,options.ascertained_position,masked_segments_dict)
        ascertained_cm1=get_cm(new_begin,new_end,begin_position,end_position,cm,recombination_rates)
        add_segment(ascertained_dict,ind_id,ascertained_cm1)
        [new_begin,new_end]=get_masked_coordinates(chromosome,options.ascertained_position,end_position,masked_segments_dict)
        ascertained_cm2=get_cm(new_begin,new_end,begin_position,end_position,cm,recombination_rates)
        add_segment(ascertained_dict,ind_id,ascertained_cm2)
        masked_sum[ind_id]+=cm-(ascertained_cm1+ascertained_cm2)
    elif IBD2=="no":
        [new_begin,new_end]=get_masked_coordinates(chromosome,begin_position,end_position,masked_segments_dict)
        new_cm=get_cm(new_begin,new_end,begin_position,end_position,cm,recombination_rates)
        add_segment(sharing_dict,ind_id,new_cm,controls)
        masked_sum[ind_id]+=cm-new_cm
    else:
        [new_begin,new_end]=get_masked_coordinates(chromosome,begin_position,end_position,masked_segments_dict)
        new_cm=get_cm(new_begin,new_end,begin_position,end_position,cm,recombination_rates)
        add_segment(ibd2_dict,ind_id,new_cm,controls)
    if controls=="yes" and options.mask_common_shared_regions!='false' and masked_segments_dict=={}:
        if chromosome not in control_segments:  
            control_segments[chromosome]=[]
        if abs(cm)>options.min_cm:
            control_segments[chromosome].append([begin_position,end_position,cm])
            # Note: this was the point at which they were being removed from unrelated founders -- default 2.5 cM 

def get_cm(begin_position,end_position,segment_begin_position,segment_end_position,cm,recombination_rates):
    global options
    if options.recombination_files is None or options.ascertained_chromosome=="no_ascertainment":
        return cm*float(end_position-begin_position)/(segment_end_position-segment_begin_position)
    else:
        for entry in recombination_rates:
            try:
                if entry.begin_position>end_position:
                    break
                if begin_position>=entry.begin_position and begin_position<=entry.end_position:
                    cm_return=entry.genetic_map_dist+entry.combined_rate*(begin_position-entry.begin_position)/1000000
                if end_position>=entry.begin_position and end_position<=entry.end_position:
                    cm_return=entry.genetic_map_dist+entry.combined_rate*(end_position-entry.begin_position)/1000000-cm_return
                    break
            except:
                msg = "%prog: Incorrect file format for recombination file"
                g.close()
                raise RuntimeError(msg)
        g.close()
        return cm_return

def get_confidence_levels(models,max_model_id,max_model_ll,confidence_statistic,model_output_file):
    global options
    global model_df
    n_0p=[]
    n_1p=[]
    n_2p=[]
    if models[max_model_id].ancestors==2 and models[max_model_id].meioses==2 and not (options.use_ibd2_siblings=="false"):
        n_2p.append(2)
    for model in models:
        if not options.model_output_file is None:
            

            if model.ancestors==2:
                dor=model.meioses-1
            else:
                dor=model.meioses
            ind_ids=model.pair_id.split(':')

            df_info = {'individual_1': ind_ids[0], 'individual_2': ind_ids[1], 'number_of_shared_ancestors': int(model.ancestors), 'degree_of_relatedness': int(dor), 'maxlnl':model.ml} 
            
            if not options.single_pair is None: # this level if-else can probably be removed now that i trim the .match file ahead of time for single_pair runs 
                
                if ind_ids[0] == options.single_pair.split(':')[0] and ind_ids[1] == options.single_pair.split(':')[1]:
                    
                    if options.write_output == True:
                        model_output_file.write(ind_ids[0]+'\t'+ind_ids[1]+'\t'+str(model.ancestors)+'\t'+str(dor)+'\t'+str(model.ml)+'\n')
                        #print (f'wrote {ind_ids[0]}-{ind_ids[1]}-{model.ancestors}-{dor} to file')

                    if options.return_output == True:
                        row_df = pd.DataFrame([df_info])
                        model_df = pd.concat([model_df, row_df], ignore_index=True)
                        #print (f'concatenated {ind_ids[0]}-{ind_ids[1]}-{model.ancestors}-{dor}')

            else:
                if options.write_output == True:
                    model_output_file.write(ind_ids[0]+'\t'+ind_ids[1]+'\t'+str(model.ancestors)+'\t'+str(dor)+'\t'+str(model.ml)+'\n')
                if options.return_output == True:
                    row_df = pd.DataFrame([df_info])
                    model_df = pd.concat([model_df, row_df], ignore_index=True)

            
        
        if model.ml+confidence_statistic>=max_model_ll:
            if model.ancestors==0:
                n_0p.append(model.meioses)
            elif model.ancestors==1:
                n_1p.append(model.meioses)
            else:
                if not (model.meioses==2 and not options.use_ibd2_siblings=="false"):
                    n_2p.append(model.meioses)
    n_0p.sort()
    n_1p.sort()
    n_2p.sort()
    if len(n_0p)>0:
        n_0p_min=n_0p[0]
        n_0p_max=n_0p[len(n_0p)-1]
    else:
        n_0p_min="none"
        n_0p_max="none"
    if len(n_1p)>0:
        n_1p_min=n_1p[0]
        n_1p_max=n_1p[len(n_1p)-1]
    else:
        n_1p_min="none"
        n_1p_max="none"
    if len(n_2p)>0:
        n_2p_min=n_2p[0]-1
        n_2p_max=n_2p[len(n_2p)-1]-1
    else:
        n_2p_min="none"
        n_2p_max="none"
    return [n_0p_min,n_0p_max,n_1p_min,n_1p_max,n_2p_min,n_2p_max]

def add_segments(beagle_marker_dict,rec_dict,chromosome_positions,sharing_dict,filename,recombination_rates,pairs,masked_segments_dict,ascertained_dict={},ibd2_dict={},control_ind=None,control_segments=None,masked_sum={}):
     global options
     '''
     reading in segment_files arg -- this is where we can find out what the .match columns are supposed to be
     '''
     
     
     ind_dict={}
     IBD2="no"
     f=open(filename,'r')

     '''handling recombination file stuff'''

     if not options.recombination_files is None:
         tmp_pos_dict={}
         for line in f.readlines():
             line_list=line.split()
             if len(line_list)!=5:
                 msg="%prog: File " + filename + " is not a Beagle fibd output file and recombination_files was specified"
                 f.close()
                 raise RuntimeError(msg)
             else:
                 chromosome=get_chromosome(filename)
                 begin_position=beagle_marker_dict[chromosome][int(line_list[2])]
                 end_position=beagle_marker_dict[chromosome][int(line_list[3])]
                 if chromosome not in chromosome_positions:
                     chromosome_positions[chromosome]={}
                 if begin_position not in chromosome_positions[chromosome]:
                     chromosome_positions[chromosome][begin_position]=0.0
                     if chromosome not in tmp_pos_dict:
                         tmp_pos_dict[chromosome]=[]
                     tmp_pos_dict[chromosome].append(begin_position)
                 if end_position not in chromosome_positions[chromosome]:
                     chromosome_positions[chromosome][end_position]=0.0
                     if chromosome not in tmp_pos_dict:
                         tmp_pos_dict[chromosome]=[]
                     tmp_pos_dict[chromosome].append(end_position)
         for chromosome,pos_list in tmp_pos_dict.iteritems():
             pos_list.sort()
             pos_index=0
             try:
                 for position in pos_list:
                     while pos_index<len(rec_dict[chromosome])-1 and rec_dict[chromosome][pos_index+1][0]<position:
                         pos_index+=1
                     if position<rec_dict[chromosome][pos_index][0]:
                         genetic_map_dist=0.0
                     else:
                         genetic_map_dist=rec_dict[chromosome][pos_index][2]+(position-rec_dict[chromosome][pos_index][0])/1e6*rec_dict[chromosome][pos_index][1]
                     chromosome_positions[chromosome][position]=genetic_map_dist
             except:
                 msg="%prog: Chromosome " + chromosome + " is not present in the recombination_files"
                 f.close()
                 raise RuntimeError(msg)

     f.close()    
     f=open(filename,'r')
     controls="no"
     if "controls" in pairs:
         if pairs["controls"]==2:
             controls="yes"
     if controls=="no":
         all_cases="no"
         if "cases" in pairs:
             if pairs["cases"]==2:
                 all_cases="yes"

     for line in f.readlines():
         line_list=line.split()

         '''handling shortened .match file'''

         if len(line_list)<12:
             if len(line_list)!=6: # changed from 5 
                 msg="%prog: File " + filename + " is not a Germline or Beagle fibd output file"
                 f.close()
                 raise RuntimeError(msg)
             else:
                '''
                there must be a working version where there's only 5 columns
                ID1, ID2, begin, end, chrom? maybe try this version 
                '''

                #ind_id=min(line_list[0],line_list[1]) + ":" + max(line_list[0],line_list[1])
                ind_id=line_list[0] + ":" + line_list[1]

                ind_dict[line_list[0]]=1
                ind_dict[line_list[1]]=1

                if controls == "yes":
                    control_ind.add(line_list[0]) 
                    control_ind.add(line_list[1]) 
                    
                #chromosome=get_chromosome(filename)
                chromosome = line_list[5] # removed the get_chromosome function since i think this version had the chr in the title instead of the file itself
                # tried updating chromosome to int instead of str
                if controls == "yes" or all_cases=="yes" or ind_id in pairs:
                    
                    '''this code that's commented out below assumes that the recombination file is passed in EVERY TIME the smaller data file is read, which isn't happening now'''
                    #begin_position=beagle_marker_dict[chromosome][int(line_list[2])]
                    #end_position=beagle_marker_dict[chromosome][int(line_list[3])]
                    #cm=chromosome_positions[chromosome][end_position]-chromosome_positions[chromosome][begin_position]

                    """we pull this data out of the GERMLINE2 output file instead of the old dicts that are otherwise empty"""                    
                    begin_position=int(line_list[2])
                    end_position=int(line_list[3])
                    cm=float(line_list[4])

                    process_segment(chromosome,ascertained_dict,sharing_dict,ibd2_dict,ind_id,cm,controls,begin_position,end_position,recombination_rates,IBD2,control_segments,masked_segments_dict,masked_sum)


         elif line_list[11]!='cM':
             msg="%prog: File " + filename + " does not provide the length of shared segments in cM"
             f.close()
             raise RuntimeError(msg)
         else:
             #print('\nProcessing .match file with >= 12 fields\n')
             '''
             handling .match file with 12 or more inputs (germline 1.5)

            0 - FID1
            1 - IID1
            2 - FID2
            3 - IID2
            4 - chrom
            5 - begin_position
            6 - end_position
            7 - rsid at beginning (not important)
            8 - rsid at end (not important)
            9 - total snps in segment (not important)
            10 - genetic length (cM)
            11 - units: cM

             '''
             ind_id=line_list[1] + ":" + line_list[3]
             ind_dict[line_list[0]]=1
             ind_dict[line_list[1]]=1
             if controls == "yes":
                 control_ind.add(line_list[1]) 
                 control_ind.add(line_list[3]) 
             if controls == "yes" or all_cases=="yes" or ind_id in pairs:
                 if len(line_list)>13 and line_list[13]=='2' and line_list[14]=='2':
                     IBD2="yes"
                 else:
                     IBD2="no"
                 cm=float(line_list[10])
                 chromosome=line_list[4]
                 begin_position=int(line_list[5])
                 end_position=int(line_list[6])
                 process_segment(chromosome,ascertained_dict,sharing_dict,ibd2_dict,ind_id,cm,controls,begin_position,end_position,recombination_rates,IBD2,control_segments,masked_segments_dict,masked_sum)
     f.close()
     for first_ind in ind_dict.keys():
         for second_ind in ind_dict.keys():
             if first_ind!=second_ind:
                 ind_id=min(first_ind,second_ind) + ":" + max(first_ind,second_ind)
                 if ind_id not in sharing_dict:
                     sharing_dict[ind_id]=[]

def shorten_match_file(pair, matchfile):
    global options

    newmatchfile = matchfile.split('.match')[0] + '_' + pair.replace(':', '-') + '.match'

    id1 = pair.split(':')[0]
    id2 = pair.split(':')[1]
    counter = 0
    with open(matchfile, 'r') as inp, open(newmatchfile, 'w+') as outp:
        for line in inp:
            if f'{id1}\t' in line and f'{id2}\t' in line:
                outp.write(line)
                counter += 1
    if options.verbose:
        print ('Length of new match file: %s lines' % str(counter))
    return newmatchfile


'''
Function that writes out a variable matrix per pair of individuals (from the pairs {} dictionary)
Matrix includes the pair of individuals and the sizes of each IBD segment shared by the pair 
'''


def return_to_primus(id1, id2, dor):

    '''
    might need to be tweaked depending on how PERL import statements work
    '''

    return (id1, id2, dor)

####################################################################################################################
####################################################################################################################
####################################################################################################################

# Wrap everything below in a function 

# we'll need to update how they call global options EVERYWHERE 

# figure out how to pull in optparse arguments from a provided dictionary 

def runner(options_arg, additional_args=None):

    global options
    argstrings = []

    if not isinstance(options_arg, Values):

        if type(options_arg) == dict:

            # instantiate default options
            parser = optparse.OptionParser()
            parser.add_option('--return_output', action="store_true",default=False, help="Returns [single pair] output in tuple format for use with PRIMUS.")      
            parser.add_option('--write_output', action="store_true",default=True, help="Writes output to .out (and/or .model) file(s).")      
            parser.add_option("--segment_files",type="string",default="*.match",help="Germline or Beagle fibd output file(s), [default: %default]")
            parser.add_option("--min_cm",type="float",default=2.5,help="minimum segment size to consider [default: %default].    If min_cm is modified, then the control_files parameter should be specified")
            parser.add_option("--max_cm",type="float",default=10.0,help="maximum segment size to consider for estimating the exponential distribution of segment sizes in the population [default: %default]")
            parser.add_option("--max_meioses",type="float",default=40,help="maximum number of meioses to consider [default: %default]")
            parser.add_option("--rec_per_meioses",type="float",default=35.2548101,help="expected number of recombination events per meioses [default: %default] from McVean et al., 2005")
            parser.add_option("--ascertained_chromosome",type="string",default="no_ascertainment",help="chromosome of ascertained disease locus")
            parser.add_option("--ascertained_position",type="int",default=-1,help="chromosomal position of ascertained disease locus")
            parser.add_option("--control_files",type="string",help="Germline or Beagle fibd output file(s) for population controls")
            parser.add_option("--control_sample_size",type="float",default=None,help="Sample size of control population.    Used only when the control_files parameter is specified, default assumes all individuals are included in the files.")
            parser.add_option("--exp_mean",type="float",default=3.197036753,help="Mean of the exponential distribution of shared segment size in the population [default: %default] from HapMap 2.0 CEU.    This parameter is ignored if mask_common_shared_regions is specified.")
            parser.add_option("--pois_mean",type="float",default=13.73,help="Mean of the Poisson distribution of the number of segments shared between a pair of individuals in the population [default: %default] from HapMap 2.0 CEU.    This parameter is ignored if mask_common_shared_regions is specified.")

            ######################
            # OLD 
            parser.add_option("--pair_file",type="string",help="Restrict pairwise comparisons to the pairs specified in this file")
            # NEW
            parser.add_option("--single_pair",type="string",help="Restrict pairwise comparisons to the pairs specified in this flag")
            ######################

            parser.add_option("--number_of_ancestors",type="int",help="Restrict relationships to [1] one parent (half-sibs/cousins), [2] two parents (full-sibs/cousins), or [0] (parent-offspring/grandparent-granchild).    Default considers all possibilities") 
            parser.add_option("--number_of_chromosomes",type="int",default=22,help="Number of chromosomes [default: %default]")
            parser.add_option("--sibling_option",type="string",default="true",help="This option was deprecated in version 1.7")
            parser.add_option("--sibling_segment_length",type="string",default="true",help="This option was deprecated in version 1.7")
            parser.add_option("--use_ibd2_siblings",type="string",default="false",help="If IBD2 data is present in the segment_file, this option will use IBD2 to detect sibling relationships.    [default: %default]")
            parser.add_option("--parent_offspring_option",type="string",default="true",help="Option to evaluate potential parent-offspring and sibling relationships based on total proportion of the genome that is shared ibd1 [default: %default]")
            parser.add_option("--parent_offspring_zscore",type="float",default=2.33,help="Zscore for rejecting a sibling relationship in favor of a parent-offspring relationship [default: %default, alpha=0.01]    Used only in combination with parent_offspring_option")
            parser.add_option("--adjust_pop_dist",type="string",default="false",help="Option to adjust the population distribution of shared segments downward for segments that could not be detected due to recent ancestry [default: %default]")
            parser.add_option("--confidence_level",type="float",default=0.95,help="Confidence level for confidence interval around the estimated degree of relationship.    If the confidence interval includes no relationship, then no_sig_rel will be reported for the estimated_degree_of_relationship [default: %default]")
            parser.add_option("--output_file",type="string",default="output/ersa.out",help="ERSA output file [default: %default]")
            parser.add_option("--mask_common_shared_regions",type="string",default="false",help="excludes chromosomal regions that are commonly shared from evaluation.    Used only when the control_files or mask_region_file parameter is specified [default: %default].")
            parser.add_option("--mask_region_cross_length",type="int",default=1000000,help="length in base pairs that a shared segment must extend past a masked segment in order to avoid truncation.    Used only when mask_common_shared_regions parameter is specified [default: %default].")
            parser.add_option("--mask_region_file",type="string",help="file containing chromosomal regions to exclude from from evaluation.    Used only when mask_common_shared_regions parameter is specified.")
            parser.add_option("--mask_region_threshold",type="float",default=4.0,help="Threshold for the ratio of observed vs. expected segment sharing in controls before a region will be masked.    Used only in conjunction with control_files and mask_common_shared_regions parameters when mask_region_file is not specified [default: %default].")
            parser.add_option("--mask_region_simulation_count",type="int",default=0,help="This option will perform simulations of the null distribution of shared segment locations in controls and will write the results of the simulations to output_file.sim.    The simulations are very slow and are not used directly in estimating relationships but allow the user to determine the max_region_threshold that meets a particular significance threshold for a given control dataset.    Used only when mask_common_shared_regions parameter is specified [default: %default].")
            parser.add_option("--recombination_files",type="string",help="file containing genetic distances for all chromosomes.    This parameter must be specified with Beagle fibd input files")
            parser.add_option("--beagle_markers_files",type="string",help="Beagle marker files (one file required for each chromosome, wildcards required, ex: chr*beagle.marker).    Each filename must begin with the chromosome name followed by a period.    This parameter must be specified with Beagle fibd input files")
            parser.add_option("--model_output_file",type="string",default=None,help="Specifies an output file to report likelihoods for all models [default: %default].")
            parser.add_option('--verbose', action="store_true", default=None, help="Determines whether or not you want to log console output print statements.")      

            options, args = parser.parse_args()

            # add dict args in here
            for key, value in options_arg.items():
                print (f'New option value read in from dict: {key}:{value}')
                argstrings.append(f'--{key}={value}')
                setattr(options, key, value)

            # print all of them
            # for option, value in vars(options).items():
            #     print(f"{option}: {value}")

        else:
            raise TypeError("Options must be in dictionary format if not provided as options via script mode. Please try again")
            sys.exit


    ################################

    global min_ll_constant
    min_ll_constant=-9999999999

    output_file=open(options.output_file,'w')

    #verbose = True if not options.verbose is None else False

    if not options.model_output_file is None:
        model_output_file=open(options.model_output_file,'w') 
    else:
        model_output_file=None

    output_file.write('# ersa version 2.2\n')

    if additional_args is not None:
        for arg in additional_args:
            output_file.write('# ' + arg + '\n')
    else:
        for arg in argstrings:
            output_file.write('# ' + arg + '\n')

    [confidence_statistic,confidence_level]=set_confidence(options.confidence_level)
    recombination_rates=[]
    rec_dict={}

    if not options.recombination_files is None:
        for recomb_filename in glob.glob(options.recombination_files):
            last=1
            if options.verbose:
                print ("Processing recombination rate file "+recomb_filename)
            g=open(recomb_filename,'r')
            line=g.readline()
            line=g.readline()
            while line:
                try:
                    line_list=line.split()
                    if line_list[0] not in rec_dict:
                        rec_dict[line_list[0]]=[]
                    rec_dict[line_list[0]].append((float(line_list[1]),float(line_list[2]),float(line_list[3])))
                    if line_list[0]==options.ascertained_chromosome:
                        recombination_rates.append(rec_entry(line_list[0],last,line_list[1],line_list[2],line_list[3]))
                    last=float(line_list[1])+1
                except:
                    msg = "%prog: Incorrect file format for recombination_file "+recomb_filename
                    g.close()
                    raise RuntimeError(msg)
                line=g.readline()
            if options.verbose:
                print ("...done")
            g.close()

    beagle_marker_dict={}
    if not options.beagle_markers_files is None:
        if options.recombination_files is None:
            msg = "%prog: beagle_marker_files parameter specified but no recombination_files provided"
            raise RuntimeError(msg)
        else:
            for bm_filename in glob.glob(options.beagle_markers_files):
                chrom_tmp=bm_filename.split(".")
                chromosome=chrom_tmp[0]
                beagle_marker_dict[chromosome]={}
                if options.verbose:
                    print ("Reading beagle marker file " + bm_filename + " for chromosome " + chromosome)
                g=open(bm_filename,'r')
                marker_index=0
                last_position=0
                for line in g.readlines():
                    try:
                        line_list=line.split()
                        beagle_marker_dict[chromosome][marker_index]=int(line_list[1])
                        last_position=int(line_list[1])
                        marker_index+=1
                    except:
                        msg="%prog: Invalid file format for beagle marker file "
                        raise RuntimeError(msg)
                beagle_marker_dict[chromosome][marker_index]=last_position+1
                if options.verbose:
                    print ("...done")


    ########################################
    # this is what we want to leverage in the library (most likely)

    if not options.pair_file is None:
        if options.verbose:
            print ("Processing pair file")
        pairs={}
        g=open(options.pair_file,'r')
        for line in g.readlines():
            try:
                line_list=line.split()
                pair_string=min(line_list[0],line_list[1])+":"+max(line_list[0],line_list[1])
                pairs[pair_string]=1
            except:
                msg = "%prog: Incorrect file format for pair_file"
                g.close()
                raise RuntimeError(msg)
        if options.verbose:
            print ("...done")
            print (pairs)


    # Single pair id1:id2 proccessing here
    else:
        if not options.single_pair is None:
            # split data
            person1 = options.single_pair.split(':')[0]
            person2 = options.single_pair.split(':')[1]
            pairstr = "%s:%s" % (person1, person2)
            pairs = {pairstr : 1}
            if options.verbose:
                print (person1, person2)
        else: # old else case
            pairs={"cases":2}

    ################
    # at this point we have a dictionary of pairs 

    ind_sharing={}
    ibd2_sharing={}
    ascertained_sharing={}
    chromosome_positions={}
    control_segments={}
    masked_segments_dict={}
    total_masked_length=0

    if options.control_files is None:
        emp_lambda=1/(options.exp_mean-options.min_cm) 
        emp_segment_lambda=options.pois_mean
        if options.mask_common_shared_regions!='false' and options.mask_region_file is None:
                msg="mask_common_shared_regions parameter must be specified with either control_files or mask_region_file parameter"
                raise RuntimeError(msg)
    else:
        cont_sharing={}
        cont_ind=set([])
        if len(glob.glob(options.control_files))==0: 
            msg="Control file " + options.control_files + " does not exist"
            raise RuntimeError(msg)
        for ctrl_filename in glob.glob(options.control_files):

            if options.verbose:
                print ("Processing control file " + ctrl_filename)

            add_segments(beagle_marker_dict,rec_dict,chromosome_positions,cont_sharing,ctrl_filename,recombination_rates,{"controls":2},masked_segments_dict,{},{},cont_ind,control_segments)
            if options.verbose:
                print ("...done")
        if options.control_sample_size is None:
            control_count=len(cont_ind)
        else:
            control_count=options.control_sample_size
        segment_counts={}
        total_segment_count=0
        total_segment_length=0.0
        emp_shared_segment_sum=0.0
        tmp_pair_count=0
        for pair_id,segments in cont_sharing.items():
            segment_count=len(segments)
            tmp_pair_count+=1
            emp_shared_segment_sum+=segment_count
            for segment in segments:
                total_segment_count+=1
                total_segment_length+=segment
        emp_segment_lambda=emp_shared_segment_sum/float(tmp_pair_count)
        output_file.write("#Mean number of shared segments between pairs of individuals in the control file(s): " + str(emp_segment_lambda) + '\n')
        exp_mean=total_segment_length/total_segment_count
        output_file.write("#Mean shared segment size in the control file(s): " + str(exp_mean) + '\n')
        emp_lambda=1/(exp_mean-options.min_cm)
        if options.mask_common_shared_regions!='false' and (options.mask_region_file is None or options.mask_region_simulation_count>0):
            if options.verbose:
                print ("Identifying common shared regions to mask")
            mask_region_segments={}
            for chromosome,segment_entries in control_segments.items():
                segment_entries.sort(key=operator.itemgetter(0))
                for i in range(len(segment_entries)):
                    expected_length=segment_entries[i][2]*(total_segment_length/(options.rec_per_meioses*100))
                    observed_length=get_total_overlap(i,segment_entries,expected_length)
                    if observed_length/expected_length>options.mask_region_threshold:
                        if chromosome not in mask_region_segments:
                            mask_region_segments[chromosome]=[]
                        mask_region_segments[chromosome].append(segment_entries[i])
            if options.verbose:
                print ("...Merging overlapping segments")
            for chromosome,masked_segments in mask_region_segments.items():
                masked_segments.sort(key=operator.itemgetter(0))
                i=0
                while i<len(masked_segments):
                    first_segment=masked_segments[i]
                    j=i+1
                    while j<len(masked_segments) and masked_segments[j][0]<=first_segment[1]:
                        first_segment[1]=max(first_segment[1],masked_segments[j][1])
                        j+=1
                    for k in range(i+1,j):
                        del masked_segments[i+1]
                    i+=1
            mask_file=open(options.output_file+'.msk','w')
            mask_file.write('chromosome\tbegin_position\tend_position\n')
            for chromosome,masked_segments in mask_region_segments.items():
                for segment in masked_segments:
                    mask_file.write(chromosome+'\t'+str(segment[0])+'\t'+str(segment[1])+'\n')
            mask_file.close()

            if options.verbose:
                print ("Masked regions written to "+options.output_file+'.msk')

            if options.mask_region_simulation_count>0:
                if options.verbose:
                    print ("Simulating null distribution of overlapping segments")
                
                segment_list=[]
                for i in range(len(segment_entries)):
                    segment_list.append([0.0,segment_entries[i][2],segment_entries[i][2],observed_length,0,chromosome,segment_entries[i][0],segment_entries[i][1],0.0])
                simulate_segments(segment_list)
                sim_file=open(options.output_file+'.sim','w')
                segment_list.sort(key=operator.itemgetter(0,6))
                sim_file.write('chromosome\tbegin_position\tend_position\tsegment_size\ttotal_observed_length_of_segment_sharing\ttotal_expected_length_of_segment_sharing\tratio_of_observed_to_expected\tnumber_of_simulations_exceeding_observed_to_expected_ratio\n')
                for segment in segment_list: 
                    expected_length=segment[2]*(total_segment_length/(options.rec_per_meioses*100))
                    sim_file.write(segment[5]+'\t'+str(segment[6])+'\t'+str(segment[7])+'\t'+str(segment[2])+'\t'+str(segment[3])+'\t'+str(expected_length)+'\t'+str(segment[3]/expected_length)+'\t'+str(segment[4])+'\n')
                sim_file.close()
                if options.verbose:
                    print ("Simulations complete.    Results written to "+options.output_file+'.sim')


    if options.mask_common_shared_regions!='false':

        if options.verbose:
            print ("Reading masked region file "+options.output_file+'.msk')

        if not options.mask_region_file is None:
            mask_file=options.mask_region_file
        else:
            mask_file=options.output_file+'.msk'
        m=open(mask_file,'r')
        masked_segments_dict={}
        line=m.readline()
        line=m.readline()
        last_segment=['chr0',0,0]
        while line and line[0]!='#':
            line_list=line.split()
            total_masked_length+=int(line_list[2])-int(line_list[1])
            if line_list[0] not in masked_segments_dict:
                masked_segments_dict[line_list[0]]=[[int(line_list[1]),int(line_list[2])]]
            elif line_list[0]==last_segment[0]:
                if int(line_list[1])<last_segment[2]:
                    msg="Mask file "+options.output_file+'.msk is not sorted or contains overlapping segments'
                    raise RuntimeError(msg)
                elif int(line_list[1])-options.mask_region_cross_length<last_segment[2]:
                    masked_segments_dict[line_list[0]][len(masked_segments_dict[line_list[0]])-1][1]=int(line_list[2])
                else:
                    masked_segments_dict[line_list[0]].append([int(line_list[1]),int(line_list[2])])
            last_segment=[line_list[0],int(line_list[1]),int(line_list[2])]
            line=m.readline()

        if options.verbose:
            print ("...done")

        if options.control_files is None:
            while line:
                line_list=line.split()
                if line_list[1]=='exp_mean':
                    if options.verbose:
                        print ("Setting exp_mean and emp_lambda paramters to value contained in    masked region file "+options.output_file+'.msk')
                    emp_lambda=1/(float(line_list[2])-options.min_cm)
                elif line_list[1]=='pois_mean':
                    emp_segment_lambda=float(line_list[2])
                line=m.readline()
            m.close()
        else:
            m.close()
            control_segments={}
            cont_sharing={}
            cont_ind=set([])
            for ctrl_filename in glob.glob(options.control_files):
                if options.verbose:
                    print ("Masking segments in control file " + ctrl_filename)
                add_segments(beagle_marker_dict,rec_dict,chromosome_positions,cont_sharing,ctrl_filename,recombination_rates,{"controls":2},masked_segments_dict,{},{},cont_ind,control_segments)
                if options.verbose:
                    print ("...done")
            if options.control_sample_size is None:
                control_count=len(cont_ind)
            else:
                control_count=options.control_sample_size
            segment_counts={}
            total_segment_count=0
            total_segment_length=0.0
            emp_shared_segment_sum=0.0
            tmp_pair_count=0
            for pair_id,segments in cont_sharing.items():
                segment_count=len(segments)
                tmp_pair_count+=1
                emp_shared_segment_sum+=segment_count
                for segment in segments:
                    total_segment_count+=1
                    total_segment_length+=segment
            mask_file=open(options.output_file+'.msk','a')
            emp_segment_lambda=emp_shared_segment_sum/float(tmp_pair_count)
            output_file.write("#Mean number of shared segments between pairs of individuals in the control file(s) after masking: " + str(emp_segment_lambda) + '\n')
            mask_file.write('# pois_mean '+str(emp_segment_lambda)+'\n')
            exp_mean=total_segment_length/total_segment_count
            mask_file.write('# exp_mean '+str(exp_mean)+'\n')
            mask_file.close()
            output_file.write("#Mean shared segment size in the control file(s) after masking: " + str(exp_mean) + '\n')
            emp_lambda=1/(exp_mean-options.min_cm)

    masked_sum={}
    if len(glob.glob(options.segment_files))==0: 
        msg="Segment file " + options.segment_files + " does not exist"
        raise RuntimeError(msg)

    for segment_filename in glob.glob(options.segment_files):
        if options.verbose:
            print ("Reading segment file " + segment_filename)

        # Restricts .match file down to only the pair we're interested in (if that flag is enabled) 
        if not options.single_pair is None:
            segment_filename = shorten_match_file(options.single_pair, segment_filename)

        add_segments(beagle_marker_dict,rec_dict,chromosome_positions,ind_sharing,segment_filename,recombination_rates,pairs,masked_segments_dict,ascertained_sharing,ibd2_sharing,set([]),{},masked_sum)
        
        if options.verbose:
            print ("...done")

    if options.verbose:
        print("Estimating recent ancestry")

    ########## Come back here after going through add_segments() function

    global genetic_map
    genetic_map=options.rec_per_meioses*100.0-(total_masked_length/1e6)*1.3

    # header

    output_file.write("individual_1\tindividual_2\test_number_of_shared_ancestors\test_degree_of_relatedness\t"+str(confidence_level)+" CI_2p_lower\t2p_upper\t1p_lower\t1p_upper\t0p_lower\t0p_upper\tmaxlnl_relatedness\tmaxlnl_unrelatedness"+'\n')

    if not options.model_output_file is None:

        if options.return_output == True:
            global model_df
            model_df = pd.DataFrame(columns=['individual_1', 'individual_2', 'number_of_shared_ancestors', 'degree_of_relatedness', 'maxlnl'])

        # instantiate df here 
        column_names = ['individual_1', 'individual_2', 'number_of_shared_ancestors', 'degree_of_relatedness', 'maxlnl']
        columnstr = ('\t'.join(column_names) + '\n')
        model_output_file.write(columnstr)

    for ind_id,ind_item in ind_sharing.items():
        print (ind_id)
        ind_item.sort(key=abs,reverse=True)
        n=1
        s=len(ind_item)
        models=[]
        while n<options.max_meioses+1:
            max_ll_0p=min_ll_constant
            max_ll_1p=min_ll_constant
            max_ll_2p=min_ll_constant
            for i in range(s+1):
                values_related=ind_item[:i]
                values_background=ind_item[i:]
                if ind_id in ascertained_sharing:
                    ascertained_values=ascertained_sharing[ind_id]
                else:
                    ascertained_values=[]
                ll=background_ll(values_background,emp_segment_lambda,emp_lambda,ascertained_values,values_related) 
                ll_0p=related_0p_ll(values_related,n)
                ll_1p=related_1p_ll(values_related,n)
                ll_2p=related_2p_ll(values_related,n)
                if len(ascertained_values)>0:
                    ll_2=background_ll(values_background,emp_segment_lambda,emp_lambda,[],values_related,ascertained_values) 
                    ll_0p_2=related_0p_ll(values_related,n,ascertained_values)
                    ll_1p_2=related_1p_ll(values_related,n,ascertained_values)
                    ll_2p_2=related_2p_ll(values_related,n,ascertained_values)
                else:
                    ll_2=min_ll_constant
                    ll_0p_2=min_ll_constant
                    ll_1p_2=min_ll_constant
                    ll_2p_2=min_ll_constant
                ll_total_0p=max(ll+ll_0p,ll_2+ll_0p_2)
                ll_total_1p=max(ll+ll_1p,ll_2+ll_1p_2)
                ll_total_2p=max(ll+ll_2p,ll_2+ll_2p_2)
                if ll_total_0p>max_ll_0p:
                    max_ll_0p=ll_total_0p
                    max_background_0p=i
                if ll_total_1p>max_ll_1p:
                    max_ll_1p=ll_total_1p
                    max_background_1p=i
                if ll_total_2p>max_ll_2p:
                    max_ll_2p=ll_total_2p
                    max_background_2p=i
            models.append(model_class(0,max_ll_0p,n,max_background_0p,s,ind_item,ind_id))
            if n!=1:
                models.append(model_class(2,max_ll_2p,n,max_background_2p,s,ind_item,ind_id))
                models.append(model_class(1,max_ll_1p,n,max_background_1p,s,ind_item,ind_id))
            n+=1
        ibd2_models=[]
        sibs_override=False
        second_relationship=False
        if options.use_ibd2_siblings=="true" and ind_id in ibd2_sharing:
            ibd1_total=0.0
            for segment in ind_item:
                ibd1_total+=segment
            ind_item_ibd2=ibd2_sharing[ind_id]
            ind_item_ibd2.sort(key=abs,reverse=True)
            n=1
            while n<options.max_meioses+1:
                max_ll_1p=min_ll_constant
                max_ll_2p=min_ll_constant
                if n==2:
                    s=len(ind_item_ibd2)
                else:
                    s=0
                unrelated_ll=background_ll(ind_item_ibd2,emp_segment_lambda*ibd1_total/genetic_map,emp_lambda)
                for i in range(s+1):
                    values_related=ind_item_ibd2[:i]
                    values_background=ind_item_ibd2[i:]
                    background_ll_ibd2=background_ll(values_background,emp_segment_lambda*ibd1_total/genetic_map,emp_lambda,[],values_related) 
                    if n==2:
                        ll_ibd2_1p=0.0
                        ll_ibd2_2p=ibd2_sib_ll(values_related)
                    elif n>=6:
                        ll_ibd2_1p=related_1p_ll(values_related,n,[],ibd1_total/genetic_map)
                        ll_ibd2_2p=0.0
                    else:
                        ll_ibd2_1p=0.0
                        ll_ibd2_2p=0.0
                    ll_total_1p=background_ll_ibd2+ll_ibd2_1p
                    ll_total_2p=background_ll_ibd2+ll_ibd2_2p
                    if ll_total_1p>max_ll_1p and i==0:
                        max_ll_1p=ll_total_1p
                        max_background_1p=i
                    if ll_total_2p>max_ll_2p:
                        max_ll_2p=ll_total_2p
                        max_background_2p=i
                ibd2_models.append(model_class(0,max_ll_1p,n,max_background_1p,s,ind_item_ibd2,ind_id))
                ibd2_models.append(model_class(1,max_ll_1p,n,max_background_1p,s,ind_item_ibd2,ind_id))
                ibd2_models.append(model_class(2,max_ll_2p,n,max_background_2p,s,ind_item_ibd2,ind_id))
                n+=1
            max_model_ll=min_ll_constant
            for model_id in range(len(ibd2_models)):
                current_model=ibd2_models[model_id]
                if options.number_of_ancestors is None or (options.number_of_ancestors==1 and current_model.ancestors==1) or (options.number_of_ancestors==2 and current_model.ancestors==2) or (options.number_of_ancestors==0 and current_model.ancestors==0):
                    if current_model.ml>max_model_ll:
                        max_model_ll=current_model.ml
                        max_model_id=model_id
                        max_model=current_model
            if max_model_ll!=min_ll_constant and max_model_ll>unrelated_ll+confidence_statistic:
                if max_model.ancestors==2 and max_model.meioses==2: 
                    sibs_override=True
                else:
                    second_relationship=max_model.meioses
        if not sibs_override:
            max_model_ll=min_ll_constant
            for model_id in range(len(models)):
                current_model=models[model_id]
                if options.number_of_ancestors is None or (options.number_of_ancestors==1 and current_model.ancestors==1) or (options.number_of_ancestors==2 and current_model.ancestors==2) or (options.number_of_ancestors==0 and current_model.ancestors==0):
                    if current_model.ml>max_model_ll and (options.use_ibd2_siblings=="false" or not (current_model.ancestors==2 and current_model.meioses==2)):
                        max_model_ll=current_model.ml
                        max_model_id=model_id
                        max_model=current_model
            if options.parent_offspring_option=="true": 
                mean=options.rec_per_meioses*100*3/4
                sd=math.sqrt((mean/100.0)*4*25**2)
                total_segment_length=0.0
                for segment in ind_item:
                    total_segment_length+=segment
                if ind_id in masked_sum:
                    total_segment_length+=masked_sum[ind_id]
                if total_segment_length>mean+options.parent_offspring_zscore*sd and (options.number_of_ancestors==0 or options.number_of_ancestors is None):
                    for model_id in range(len(models)):
                        if models[model_id].ancestors==0 and models[model_id].meioses==1:
                            max_model_id=model_id
                            max_model_ll=models[model_id].ml
                            max_model=models[model_id]
                elif options.use_ibd2_siblings=="false" and total_segment_length/(options.rec_per_meioses*100)>0.625 and (options.number_of_ancestors==2 or options.number_of_ancestors is None):
                    for model_id in range(len(models)):
                        if models[model_id].ancestors==2 and models[model_id].meioses==2:
                            max_model_id=model_id
                            max_model_ll=models[model_id].ml
                            max_model=models[model_id]
        if sibs_override: 
            [n_0p_min,n_0p_max,n_1p_min,n_1p_max,n_2p_min,n_2p_max]=get_confidence_levels(ibd2_models,max_model_id,max_model_ll,confidence_statistic,model_output_file)
            total_expected_ibd1=min((emp_segment_lambda*emp_lambda)/(genetic_map),1.0)
            ll=background_ll(values_background,emp_segment_lambda*total_expected_ibd1,emp_lambda)
        else: 
            [n_0p_min,n_0p_max,n_1p_min,n_1p_max,n_2p_min,n_2p_max]=get_confidence_levels(models,max_model_id,max_model_ll,confidence_statistic,model_output_file)
            ll=background_ll(ind_item,emp_segment_lambda,emp_lambda,ascertained_values)
        [ind1,ind2]=ind_id.split(":")
        

        # right here we should be able to check ind1 and ind2 if they match the one in pairs dict 
        
        if not options.single_pair is None:
            if ind1 == person1 and ind2 == person2: # finding the correct pair to write out

                # no significant relatedness
                if ll+confidence_statistic>=max_model_ll:
                    dor = 'no_sig_rel'
                    output_file.write(ind1+"\t" + ind2+"\t0\tno_sig_rel\t"+str(n_2p_min)+"\t"+str(n_2p_max)+"\t"+str(n_1p_min)+"\t"+str(n_1p_max)+"\t"+str(n_0p_min)+'\t'+str(n_0p_max)+'\t'+str(max_model_ll)+"\t"+str(ll)+"\n")
                
                # there is some relatedness
                else:
                    if max_model.ancestors==2:
                        dor=str(max_model.meioses-1) # dor = degree of relatedness
                    else:
                        dor=str(max_model.meioses)
                    if second_relationship:
                        dor+="("+second_relationship+")"
                    output_file.write(ind1+"\t"+ind2+"\t"+str(int(max_model.ancestors))+"\t"+dor+"\t"+str(n_2p_min)+"\t"+str(n_2p_max)+"\t"+str(n_1p_min)+"\t"+str(n_1p_max)+"\t"+str(n_0p_min)+'\t'+str(n_0p_max)+'\t'+str(max_model_ll)+"\t"+str(ll)+"\n")


        
        else: # Writing every pair, not a single one 

            if ll+confidence_statistic>=max_model_ll:
                output_file.write(ind1+"\t" + ind2+"\t0\tno_sig_rel\t"+str(n_2p_min)+"\t"+str(n_2p_max)+"\t"+str(n_1p_min)+"\t"+str(n_1p_max)+"\t"+str(n_0p_min)+'\t'+str(n_0p_max)+'\t'+str(max_model_ll)+"\t"+str(ll)+"\n")
            else:
                if max_model.ancestors==2:
                    dor=str(max_model.meioses-1)
                else:
                    dor=str(max_model.meioses)
                if second_relationship:
                    dor+="("+second_relationship+")"
                output_file.write(ind1+"\t"+ind2+"\t"+str(int(max_model.ancestors))+"\t"+dor+"\t"+str(n_2p_min)+"\t"+str(n_2p_max)+"\t"+str(n_1p_min)+"\t"+str(n_1p_max)+"\t"+str(n_0p_min)+'\t'+str(n_0p_max)+'\t'+str(max_model_ll)+"\t"+str(ll)+"\n")

    '''Add in support for building a dataframe instead of file writing if the argument is passed'''


    if options.verbose:
        print ("ERSA completed successfully\n")

    # running ersa pairwise generates a filtered .match file for the pair, which is unnecessary, so i'm auto-deleting them after running
    if not options.single_pair is None:
        os.system(f'rm {segment_filename}')

    if options.return_output == True:
        return model_df


if __name__ == "__main__":

    # run-time parameters ###################################################################################################
    parser=optparse.OptionParser()

    parser.add_option('--return_output', action="store_true",default=False, help="Returns [single pair] output in tuple format for use with PRIMUS.")      
    parser.add_option('--write_output', action="store_true",default=True, help="Writes output to .out (and/or .model) file(s).")      
    parser.add_option("--segment_files",type="string",default="*.match",help="Germline or Beagle fibd output file(s), [default: %default]")
    parser.add_option("--min_cm",type="float",default=2.5,help="minimum segment size to consider [default: %default].    If min_cm is modified, then the control_files parameter should be specified")
    parser.add_option("--max_cm",type="float",default=10.0,help="maximum segment size to consider for estimating the exponential distribution of segment sizes in the population [default: %default]")
    parser.add_option("--max_meioses",type="float",default=40,help="maximum number of meioses to consider [default: %default]")
    parser.add_option("--rec_per_meioses",type="float",default=35.2548101,help="expected number of recombination events per meioses [default: %default] from McVean et al., 2005")
    parser.add_option("--ascertained_chromosome",type="string",default="no_ascertainment",help="chromosome of ascertained disease locus")
    parser.add_option("--ascertained_position",type="int",default=-1,help="chromosomal position of ascertained disease locus")
    parser.add_option("--control_files",type="string",help="Germline or Beagle fibd output file(s) for population controls")
    parser.add_option("--control_sample_size",type="float",default=None,help="Sample size of control population.    Used only when the control_files parameter is specified, default assumes all individuals are included in the files.")
    parser.add_option("--exp_mean",type="float",default=3.197036753,help="Mean of the exponential distribution of shared segment size in the population [default: %default] from HapMap 2.0 CEU.    This parameter is ignored if mask_common_shared_regions is specified.")
    parser.add_option("--pois_mean",type="float",default=13.73,help="Mean of the Poisson distribution of the number of segments shared between a pair of individuals in the population [default: %default] from HapMap 2.0 CEU.    This parameter is ignored if mask_common_shared_regions is specified.")

    ######################
    # OLD 
    parser.add_option("--pair_file",type="string",help="Restrict pairwise comparisons to the pairs specified in this file")
    # NEW
    parser.add_option("--single_pair",type="string",help="Restrict pairwise comparisons to the pairs specified in this flag")
    ######################

    parser.add_option("--number_of_ancestors",type="int",help="Restrict relationships to [1] one parent (half-sibs/cousins), [2] two parents (full-sibs/cousins), or [0] (parent-offspring/grandparent-granchild).    Default considers all possibilities") 
    parser.add_option("--number_of_chromosomes",type="int",default=22,help="Number of chromosomes [default: %default]")
    parser.add_option("--sibling_option",type="string",default="true",help="This option was deprecated in version 1.7")
    parser.add_option("--sibling_segment_length",type="string",default="true",help="This option was deprecated in version 1.7")
    parser.add_option("--use_ibd2_siblings",type="string",default="false",help="If IBD2 data is present in the segment_file, this option will use IBD2 to detect sibling relationships.    [default: %default]")
    parser.add_option("--parent_offspring_option",type="string",default="true",help="Option to evaluate potential parent-offspring and sibling relationships based on total proportion of the genome that is shared ibd1 [default: %default]")
    parser.add_option("--parent_offspring_zscore",type="float",default=2.33,help="Zscore for rejecting a sibling relationship in favor of a parent-offspring relationship [default: %default, alpha=0.01]    Used only in combination with parent_offspring_option")
    parser.add_option("--adjust_pop_dist",type="string",default="false",help="Option to adjust the population distribution of shared segments downward for segments that could not be detected due to recent ancestry [default: %default]")
    parser.add_option("--confidence_level",type="float",default=0.95,help="Confidence level for confidence interval around the estimated degree of relationship.    If the confidence interval includes no relationship, then no_sig_rel will be reported for the estimated_degree_of_relationship [default: %default]")
    parser.add_option("--output_file",type="string",default="output/ersa.out",help="ERSA output file [default: %default]")
    parser.add_option("--mask_common_shared_regions",type="string",default="false",help="excludes chromosomal regions that are commonly shared from evaluation.    Used only when the control_files or mask_region_file parameter is specified [default: %default].")
    parser.add_option("--mask_region_cross_length",type="int",default=1000000,help="length in base pairs that a shared segment must extend past a masked segment in order to avoid truncation.    Used only when mask_common_shared_regions parameter is specified [default: %default].")
    parser.add_option("--mask_region_file",type="string",help="file containing chromosomal regions to exclude from from evaluation.    Used only when mask_common_shared_regions parameter is specified.")
    parser.add_option("--mask_region_threshold",type="float",default=4.0,help="Threshold for the ratio of observed vs. expected segment sharing in controls before a region will be masked.    Used only in conjunction with control_files and mask_common_shared_regions parameters when mask_region_file is not specified [default: %default].")
    parser.add_option("--mask_region_simulation_count",type="int",default=0,help="This option will perform simulations of the null distribution of shared segment locations in controls and will write the results of the simulations to output_file.sim.    The simulations are very slow and are not used directly in estimating relationships but allow the user to determine the max_region_threshold that meets a particular significance threshold for a given control dataset.    Used only when mask_common_shared_regions parameter is specified [default: %default].")
    parser.add_option("--recombination_files",type="string",help="file containing genetic distances for all chromosomes.    This parameter must be specified with Beagle fibd input files")
    parser.add_option("--beagle_markers_files",type="string",help="Beagle marker files (one file required for each chromosome, wildcards required, ex: chr*beagle.marker).    Each filename must begin with the chromosome name followed by a period.    This parameter must be specified with Beagle fibd input files")
    parser.add_option("--model_output_file",type="string",default=None,help="Specifies an output file to report likelihoods for all models [default: %default].")
    parser.add_option('--verbose', action="store_true", default=False, help="Determines whether or not you want to log console output print statements.")      

    (options, sys.args) = parser.parse_args()

    runner(options, sys.args)

    # If it makes it this far, it means that the options are being passed to the runner in optparse format