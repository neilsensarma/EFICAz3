import os
import subprocess
from multiprocessing import Pool
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices
import shutil

'''
Use when working locally

input_path = '/Users/neil/Desktop/EFICAz3/eficaz3/data/DATASET/UNIPROT_SPROT/4EC/'
output_path = '/Users/neil/Desktop/EFICAz3/eficaz3/data/DATASET/UNIPROT_SPROT/4EC_HMM_trial2/'
'''

'''
Use when working on the server
'''
input_path = '/storage/home/hhive1/nsensarma3/data/eficaz3/data/DATASET/UNIPROT_SPROT/3EC/'
output_path = '/storage/home/hhive1/nsensarma3/data/eficaz3/data/DATASET/UNIPROT_SPROT/3EC_HMM_trial/'


'''
Use when working locally on 3EC files
input_path = '/Users/neil/Desktop/EFICAz3/eficaz3/data/DATASET/UNIPROT_SPROT/3EC/'
output_path = '/Users/neil/Desktop/EFICAz3/eficaz3/data/DATASET/UNIPROT_SPROT/3EC_HMM_trial/'
'''

ecs = ['uniprot_sprot_put_2.7.12','uniprot_sprot_put_1.1.99', 'uniprot_sprot_put_5.3.4']
fasta_files = [f for f in os.listdir(input_path) if f in ecs]

if os.path.exists(output_path):
    cmd = f'rm -rf {output_path}'
    os.system(cmd)

if not os.path.exists(output_path):
    os.mkdir(f'{output_path}')

aligner = PairwiseAligner(mode='global')
matrix = substitution_matrices.load('BLOSUM62')
aligner.substitution_matrix = matrix


def chiefc(fasta):
    '''
    Description: main driver code to perform enzyme family classification
    Input:
        fasta (file) = EC fasta file
    Output:
        Enzyme families (files)
    '''
    #get ec_number
    ec = fasta.split('_')[-1]
    
    #read file and get master sequence list
    master_sequences = read(fasta)
    
    #Cluster the current fasta file
    cluster_file = cluster_seqs(ec, fasta)

    #parse cluster file
    master_clusters, clusters = parse_clusters(master_sequences, cluster_file)

    #get the list of clusters sorted by sequence identity in descending order
    cluster_seq_ids = calc_seq_id(clusters)
    master_log_file = f'{output_path}{ec}.log'
    master_log = open(os.path.join(output_path, master_log_file), 'w+')
    total_sequences = master_sequences.copy() #list to reduce when sequences get hit as significant.
    i = 0 #counter to loop through the cluster_seq_ids
    all_seqs_classified = False #flag to check if all sequences have been classified
    while not all_seqs_classified:
        if total_sequences and i<=len(cluster_seq_ids)-1:
            #choose the ith cluster
            cluster = cluster_seq_ids[i]
            master_log.write(f'\nChosen cluster: {cluster}, {i+1}th iteration\n')
            
            #get sequences of the ith cluster
            subgroup_sequences = get_cluster_sequences(master_clusters, cluster)
            
            #filter sequences based on a threshold
            filtered_sequences = filter_sequences(subgroup_sequences)
            master_log.write(f'{len(subgroup_sequences)}, {len(filtered_sequences)}\n')
            total_sequences = reduce_total_sequences(total_sequences, filtered_sequences, master_sequences)
            master_log.write(f'Remaining Sequences: {len(total_sequences)}\n')
            
            #define file names for MSA and HMM and output family name
            seqs_for_msa = f'{output_path}{ec}_subgroup_for_msa{i+1}'
            msa = f'{output_path}{ec}.nr.msa{i+1}' #output MSA file 
            hmm = f'{output_path}{ec}.nr.hmm{i+1}' #output HMM file
            seqs_for_hmmsearch = f'{output_path}{ec}_rem_sequences{i+1}' #sequences to perform hmmsearch on
            fam_name = f'{output_path}{ec}.fam{i+1}' #output family name
            
            significant_hits = True #flag to check if HMMSearch returns any significant hits or not
            while significant_hits:
                #perform MSA on filtered sequences
                write_sequences(seqs_for_msa, filtered_sequences)
                msa_success = make_msa(seqs_for_msa, msa)
                
                master_log.write(f'{msa_success}\n')
                #create HMM using the MSA
                hmm_success = make_hmm(hmm, msa)
                master_log.write(f'{hmm_success}\n')
                
                #HMMSearch
                write_sequences(seqs_for_hmmsearch, total_sequences)
                hmmsearch_res, hmmsearch_success = perform_hmmsearch(ec, hmm, seqs_for_hmmsearch, i) #perform HMMSearch
                master_log.write(f'{hmmsearch_success}\n')
                
                #Parse HMMSearch results
                hits, hits_file = parse_hmmsearch(ec, hmmsearch_res, master_sequences, i)
                if hits:
                    write_sequences(hits_file, hits)
                    filtered_sequences = concatenate_seqs(ec, seqs_for_msa, hits_file, i)
                    total_sequences = reduce_total_sequences(total_sequences, hits, master_sequences)
                    master_log.write(f'No. of sequences remaning after current iteration: {len(total_sequences)}\n')
                    if not total_sequences:
                        master_log.write('All Sequences Significant!\n')
                        write_sequences(fam_name, filtered_sequences)
                        significant_hits = False
                else:
                    write_sequences(fam_name, filtered_sequences)
                    significant_hits = False #no significant hits found
                    i+=1 #move on to the next cluster
                    master_log.write(f'Number of sequences left: {len(total_sequences)}\n')
        elif total_sequences and i >= len(cluster_seq_ids):
            num = len(total_sequences)
            for j in range(num):
                fam_num = i+j+1
                file_name = f'{output_path}{ec}.fam{fam_num}'
                write_sequences(file_name, [total_sequences[j]])
            total_sequences.clear()
        else:
            master_log.write('All Sequences Classified!\n')
            all_seqs_classified = True #if all sequences have been classified, then end the loop
    cleanup(ec)
    master_log.close() #close the master log file

def read(filename):
    '''
    Description: Read the EC sequence file and retrieve the sequences and sequence ID's
    Input: 
        filename (file) = EC sequence file
    Output: 
        sequences (list) = Master list of Sequence ID's and Sequences
    '''
    sequences = [record for record in SeqIO.parse(f'{input_path}{filename}', format='fasta')]
    return sequences

def run_commands(command, log_file=None):
    if log_file:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        stdout, stderr = process.communicate()
        with open(log_file, 'w+') as log:
            log.write('Standard Output:\n')
            log.write(stdout)
            log.write('Standard Erorr:\n')
            log.write(stderr)
    else:
        subprocess.run(command, shell=True)

def cluster_seqs(ec, filename):
    '''
    Description: cluster the EC file
    Input:
        ec (int) = EC number
        filename (file) = EC sequence file
    Output:
        cluster_file (file) = file containing clusters of sequences
    '''
    cluster_file = f'{ec}.nr'
    cmd = f'cd-hit -i {input_path}{filename} -o {output_path}{cluster_file} -c 0.4 -n 2'
    cluster_log_file = f'{output_path}{ec}_cluster.log'
    run_commands(cmd, cluster_log_file)
    return cluster_file

def parse_clusters(master_sequences, cluster_file):
    '''
    Description: Parse the cluster file, for each cluster, retrieve the sequence ID's and the sequences
    Input:
        master_sequences (list) = master sequences to retrieve the sequences as cluster_file has only IDs
        cluster_file (file) = cluster file
    Output:
        master_clusters (dict) = dictionary of all clusters and their sequences and IDs
        clusters (dict) = dictionary of clusters and sequences
    '''
    master_clusters, clusters = {}, {}
    new_sequences = master_sequences.copy()
    with open(f'{output_path}{cluster_file}.clstr', 'r+') as cluster_handle:
        #print('Read the cluster file')
        lines = cluster_handle.read().strip().split('\n')  
    for i in lines:
        if i.startswith('>'):
            #print('Cluster found')
            cluster_num = f'Cluster {i.split(" ")[1]}'
            if cluster_num not in clusters:
                clusters[cluster_num], master_clusters[cluster_num] = [], []
                #print(clusters)
        elif not i.startswith('>'):            
            seq_id = i.split(' ')[1].replace('.', '').replace('>', '')
            #print(seq_id)
            for j in new_sequences:
                if j.id == seq_id:
                    temp = {f'{seq_id}': j.seq}
                    master_clusters[cluster_num].append(temp)
                    clusters[cluster_num].append(j.seq)
                    new_sequences = [record for record in new_sequences if record.id != j.id]
    return master_clusters, clusters

def calc_seq_id(clusters):
    '''
    Description: calculate the similarity ratio for each cluster
    Input: 
        clusters (dict) = dictionary of clusters and sequences
    Output:
        sorted_keys (list) = list of clusters in descending order of 
        their sequence similarity scores
    '''
    cluster_seq_ids = {}
    for k,v in clusters.items():
        num_sequences = len(v)
        total_identity = 0.0
        if num_sequences > 1:
            for i in range(num_sequences-1):
                for j in range(i+1, num_sequences):
                    seq_a, seq_b = v[i], v[j]
                    alignment = aligner.align(seq_a, seq_b)
                    similarity = alignment.score/(max(len(seq_a), len(seq_b)))
                    total_identity += similarity
        elif num_sequences <= 1:
            continue
        average_identity = total_identity/(num_sequences * (num_sequences-1)/2)
        ratio = num_sequences/average_identity
        cluster_seq_ids[k] = ratio
    sorted_clusters = dict(sorted(cluster_seq_ids.items(), key=lambda kv: kv[1], reverse=True))
    return list(sorted_clusters.keys())

def get_cluster_sequences(master_cluster_list, cluster_num):
    '''
    Description: get the sequences for each cluster
    Input:
        master_cluster_list (dict) = dictionary that contains the sequences and IDs for current cluster
        cluster_num (int) = current chosen cluster
    Output:
        sequence_records (list) = sequences for each sequence in the current cluster
    '''
    sequence_records = []
    if cluster_num in master_cluster_list:
        for record in master_cluster_list[cluster_num]:
            for seq_id, sequence in record.items():
                records = SeqRecord(id=seq_id, seq=sequence, description='')
                sequence_records.append(records)
    return sequence_records

def filter_sequences(sequence_records):
    '''
    Description: filter sequences based on threshold (85%)
    Input:
        sequence_records (list) = sequences of current cluster
    Output:
        filtered_sequences (list) = sequences that share less than 85% similarity.
    '''
    num_sequences = len(sequence_records)
    threshold = 0.85
    filtered_sequences = []
    if num_sequences <=2:
        filtered_sequences = sequence_records
        return filtered_sequences
    for i in range(num_sequences-1):
        keep_sequence = True
        for j in range(i+1, num_sequences):
            seq_a, seq_b = sequence_records[i].seq, sequence_records[j].seq
            alignment = aligner.align(seq_a, seq_b)
            similarity = (alignment.score/max(len(seq_a), len(seq_b)))/100
            if similarity >= threshold:
                keep_sequence = False
                break
        if keep_sequence:
            filtered_sequences.append(sequence_records[i])
    return filtered_sequences

def write_sequences(file_name, sequence_records):
    '''
    Description: write the provided sequences into a file
    Input: 
        file_name (str) = name of the file (includes output path)
        sequence_records (list) = SeqRecord objects
    Output:
        writes the provided sequences to the file
    '''
    with open(f'{file_name}', 'w+') as file:
        SeqIO.write(sequence_records, file, format='fasta')

def make_msa(sequences, msa_file_name):
    '''
    Description: perform MSA on the provided sequence file
    Input:
        sequences (file) = fasta file that contains sequences to perform MSA
        msa_file_name (str) = output MSA file name in FASTA format
    Output:
        performs MSA and returns MSA file.
    '''
    msa_command = f'clustalw2 -infile={sequences} -outfile={msa_file_name} -output="FASTA" -type="PROTEIN" -align -quiet -pwmatrix="BLOSUM"'
    msa_log_file = f'{msa_file_name}.log'
    #os.system(f'{msa_command} |& tee {msa_log_file}')
    run_commands(msa_command, msa_log_file)
    if os.path.isfile(msa_file_name):
        return 'MSA completed'
    
def make_hmm(hmm_file_name, msa):
    '''
    Description: build HMM using the provided MSA
    Input: 
        hmm_file_name (str) = name of the output HMM file
        msa (file) = MSA file in FASTA format
    Output:
        builds the HMM.
    '''
    hmm_command = f'hmmbuild {hmm_file_name} {msa}'
    hmm_log_file = f'{hmm_file_name}.log'
    #os.system(f'{hmm_command} |& tee {hmm_log_file}')
    run_commands(hmm_command, hmm_log_file)
    if os.path.isfile(hmm_file_name):
        return 'HMM Built'
    
def perform_hmmsearch(ec, hmm, seqs_file, cluster_num):
    '''
    Description: perform hmmsearch of a query HMM (make_hmm output) and the target database
    Input:
        ec (str) = EC number 
        hmm (file) = query HMM file
        seqs_file (file) = target database
        cluster_num (int) = current cluster
    Output:
        output (str) = results of the hmmsearch
    '''
    output = f'{output_path}{ec}.hmmsearch{cluster_num+1}'
    hmmsearch_command = f'hmmsearch {hmm} {seqs_file} > {output}'
    #os.system(hmmsearch_command)
    run_commands(hmmsearch_command)
    if os.path.isfile(f'{output}'):
        success = 'HMMSearch Completed'
    return output, success

def concatenate_seqs(ec, file1, file2, cluster_num):
    '''
    Description: combine the sequence files
    Input: 
        ec (str) = EC number
        file1 (file) = first file to combine (FASTA format)
        file2 (file) = second file to combine (FASTA format)
        cluster_num (int) = current cluster
    Output:
        (list) =  SeqRecords of the concatenated FASTA file
    '''
    combined = f'{output_path}{ec}_combined{cluster_num+1}'
    cmd = f'cat {file1} {file2} > {combined}'
    run_commands(cmd)
    return [record for record in SeqIO.parse(f'{combined}', format='fasta')]

def parse_hmmsearch(ec, hmm_res, master_sequences, cluster_num):
    '''
    Description: parse the HMMSearch results file to search for significant hits
    Input:
        ec (str) = EC number
        hmm_res (file) = resulting hmmsearch file
        master_sequences (list) = SeqRecord objects of all sequences in the EC group.
        cluster_num (int) = current cluster
    Output:
        significant_seqs (list/False) = SeqRecord objects of significant hits/False in case of no hits
        hits_file (file) = output hits file to write the significant hits as FASTA in future
    '''
    hits, significant_seqs = {}, []
    hits_file = f'{output_path}{ec}.hits{cluster_num+1}'
    with open(f'{hmm_res}', 'r+') as file:
        lines = file.read().strip().split('\n')
    for i in range(14, len(lines)):
        if lines[i].startswith('Domain'):
            break
        e_val, id = lines[i].strip().split(' ')[0], lines[i].strip().split(' ')[-1]
        if e_val == '[No': #this checks for the presence of no significant hits
            return False, False
        if id:
            hits[id] = float(e_val)
    for id, e_val in hits.items():
        if e_val < 0.05:
            for j in master_sequences:
                if id == j.id:
                    seq_rec = SeqRecord(id = j.id, seq=j.seq, description = '')
                    significant_seqs.append(seq_rec)
    return significant_seqs, hits_file

def reduce_total_sequences(ts, seqs_to_remove, master_sequences):
    '''
    Description: reduce the list of sequences that is the target database after each iteration.
    Input:
        ts (list) = current total sequences 
        seqs_to_remove (list) = sequences to remove from the total_sequences
        master_sequences (list) = master list of all sequences and their IDs in the EC group
    Output:
        ts_seqs (list) = reduced set of sequences
    '''
    t_s = [f.id for f in ts]
    h_s = [f.id for f in seqs_to_remove]
    r_s = [f for f in t_s if f not in h_s]
    ts_seqs = []
    for id in r_s:
        for record in master_sequences:
            if id == record.id:
                seq_rec = SeqRecord(id = record.id, seq=record.seq, description='')
                ts_seqs.append(seq_rec)
    return ts_seqs

def file_move(file, path=False):
    if os.path.isfile(file):
        if path:
            shutil.move(file, path)
    if os.path.isdir(file):
        pass
        
def cleanup(ec):
    '''
    Description: cleanup the excessive files being generated during this iterative procedure
    '''
    path = os.path.join(output_path, ec)
    os.makedirs(path)
    files_to_move = [f for f in os.listdir(output_path)]
    for f in files_to_move:
        if f'{ec}' in f:
            file_move(os.path.join(output_path, f), path)
        if '.dnd' in f:
            os.remove(os.path.join(output_path, f))

for fasta in fasta_files:
    chiefc(fasta)
    