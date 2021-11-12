#! /usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
from scipy.optimize import newton
import argparse
import os
import shutil
import fnmatch
import sys
import errno
import pandas as pd
import subprocess
from subprocess import call, check_output, STDOUT
import multiprocessing as mp

__version__ = 'skmer 3.2.1'

# Hard-coded param
coverage_threshold = 5
error_rate_threshold = 0.03
seq_len_threshold = 2000
default_error_rate = 0.01


def sequence_stat(sequence):
    total_length = 0
    n_reads = 0
    max_length = 0
    comp_stdout = check_output(["seqtk", "comp", sequence], stderr=STDOUT, universal_newlines=True)
    reads_stat = comp_stdout.split('\n')
    for stat in reads_stat:
        if not stat.strip():
            continue
        read_length = sum([int(x) for x in stat.split('\t')[2:6]])
        total_length += read_length
        max_length = max(max_length, read_length)
        n_reads += 1
    return int(round(1.0 * total_length / n_reads)), max_length, total_length, n_reads


def sample_reads(sequence, seed, bl_sz, bs_dir):
   
    try:
        os.makedirs(bs_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    bs_rep = os.path.join(bs_dir, os.path.split(sequence)[-1])
  
    with open(bs_rep, 'w') as fp: 
        subprocess.run(["seqtk", "sample",  "-s", str(seed), sequence, str(bl_sz)], stdout=fp) 

    return 


def cov_temp_func(x, r, p, k, l):
    lam = x * (1.0 * (l - k)) / l
    return lam * (p ** 2) * np.exp(-lam * p) - 2 * r * (p * np.exp(-lam * p) + 1 - p)


def estimate_cov(sequence, lib, k, e, nth):
    sample = os.path.basename(sequence).rsplit('.f', 1)[0]
    sample_dir = os.path.join(lib, sample)
    try:
        os.makedirs(sample_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise
    info_file = os.path.join(sample_dir, sample + '.dat')

    (l, ml, tl, n_reads) = sequence_stat(sequence)
    if ml > seq_len_threshold:
        cov = "NA"
        g_len = tl
        eps = 0
        l = "NA"
        with open(info_file, mode='w') as f:
            f.write('coverage\t{0}\n'.format(cov) + 'genome_length\t{0}\n'.format(g_len) +
                    'error_rate\t{0}\n'.format(eps) + 'read_length\t{0}\n'.format(l))
        return sample, cov, g_len, eps, l

    mercnt = os.path.join(sample_dir, sample + '.jf')
    histo_file = os.path.join(sample_dir, sample + '.hist')
    call(["jellyfish", "count", "-m", str(k), "-s", "100M", "-t", str(nth), "-C", "-o", mercnt, sequence],
         stderr=open(os.devnull, 'w'))
    histo_stderr = check_output(["jellyfish", "histo", "-h", "1000000", mercnt], stderr=STDOUT, universal_newlines=True)
    with open(histo_file, mode='w') as f:
        f.write(histo_stderr)
    os.remove(mercnt)
    count = [0]
    ksum = 0
    for item in histo_stderr.split('\n')[:-1]:
        count.append(int(item.split()[1]))
        ksum += int(item.split()[0]) * int(item.split()[1])
    if len(count) < 3:
        sys.stderr.write('Coverage of {0} is too low, not able to estimate it; no correction applied\n'.format(sample))
        cov = "NA"
        g_len = "NA"
        eps = "NA"
        with open(info_file, mode='w') as f:
            f.write('coverage\t{0}\n'.format(cov) + 'genome_length\t{0}\n'.format(g_len) +
                    'error_rate\t{0}\n'.format(eps) + 'read_length\t{0}\n'.format(l))
        return sample, cov, g_len, eps, l
    ind = min(count.index(max(count[2:])), len(count) - 2)
    if e is not None:
        eps = e
        p0 = np.exp(-k * eps)
        if ind < 2:
            r21 = 1.0 * count[2] / count[1]
            cov = newton(cov_temp_func, 0.05, args=(r21, p0, k, l))
        else:
            cov = (1.0 / p0) * (1.0 * l / (l - k)) * (ind + 1) * count[ind + 1] / count[ind]
    elif ind < 2:
        sys.stderr.write('Not enough information to co-estimate coverage and error rate of {0}; '.format(sample) +
                         'Using default error rate {0}\n'.format(default_error_rate))
        eps = default_error_rate
        p0 = np.exp(-k * eps)
        r21 = 1.0 * count[2] / count[1]
        cov = newton(cov_temp_func, 0.05, args=(r21, p0, k, l))
    else:
        gam = 1.0 * (ind + 1) * count[ind + 1] / count[ind]
        lam = (np.exp(-gam) * (gam ** ind) / np.math.factorial(ind)) * count[1] / count[ind] + gam * (1 - np.exp(-gam))
        eps = 1 - (gam / lam) ** (1.0 / k)
        cov = (1.0 * l / (l - k)) * lam
    tot_seq = 1.0 * ksum * l / (l - k)
    g_len = int(tot_seq / cov)

    if eps > error_rate_threshold or eps < 0:
        cov = "NA"
        g_len = "NA"
        eps = "NA"
        with open(info_file, mode='w') as f:
            f.write('coverage\t{0}\n'.format(cov) + 'genome_length\t{0}\n'.format(g_len) +
                    'error_rate\t{0}\n'.format(eps) + 'read_length\t{0}\n'.format(l))
        return sample, cov, g_len, eps, l

    with open(info_file, mode='w') as f:
        f.write('coverage\t{0}\n'.format(repr(cov)) + 'genome_length\t{0}\n'.format(g_len) +
                'error_rate\t{0}\n'.format(repr(eps)) + 'read_length\t{0}\n'.format(l))
    return sample, cov, g_len, eps, l


def estimate_stats(sequence, nth):
    sample = os.path.basename(sequence).rsplit('.f', 1)[0]

    (l, ml, tl, n_reads) = sequence_stat(sequence)
    if ml > seq_len_threshold:
        cov = "NA"
        g_len = tl
        eps = 0
        l = "NA"
    else:
       # Set to dummy values for reads to initialize dictionaries. 
       # Will be recomputed for each subsample.
        cov = 0.0
        g_len = tl
        eps = 0.0
        l = l
    return sample, cov, g_len, eps, l, n_reads

def create_sketch_dir(sequence, lib, ce, ge, ee, le,  nth):
    sample = os.path.basename(sequence).rsplit('.f', 1)[0]
    sample_dir = os.path.join(lib, sample)
    try:
        os.makedirs(sample_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise
    info_file = os.path.join(sample_dir, sample + '.dat')

    cov = ce[sample]
    g_len = ge[sample]
    eps = ee[sample]
    l = le[sample]
    with open(info_file, mode='w') as f:
        f.write('coverage\t{0}\n'.format(cov) + 'genome_length\t{0}\n'.format(g_len) +
                'error_rate\t{0}\n'.format(eps) + 'read_length\t{0}\n'.format(l))
    return


def sketch(sequence, lib, ce, ee, k, s, cov_thres, seed):
    sample = os.path.basename(sequence).rsplit('.f', 1)[0]
    sample_dir = os.path.join(lib, sample)
    msh = os.path.join(sample_dir, sample)
    cov = ce[sample]
    eps = ee[sample]
    if cov == "NA" and eps == 0:
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-S", str(seed), "-o", msh, sequence], stderr=open(
            os.devnull, 'w'))
        return
    elif eps == "NA":
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-S", str(seed), "-r", "-o", msh, sequence], stderr=open(
            os.devnull, 'w'))
        return
    copy_thres = int(cov / cov_thres) + 1
    if cov < cov_thres or eps == 0.0:
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-S", str(seed), "-r", "-o", msh, sequence], stderr=open(
            os.devnull, 'w'))
    else:
        call(["mash", "sketch", "-m", str(copy_thres), "-k", str(k), "-s", str(s), "-S", str(seed), "-o", msh,
              sequence], stderr=open(os.devnull, 'w'))
    return


def jacc2dist(j, k, gl1, gl2, len_penalty):
    if len_penalty:
        return 1 - (1.0 * (gl1 + gl2) * j / (1.0 * (gl1 + gl2) * (1 + j) / 2)) ** (1.0 / k)
    else:
        return 1 - (1.0 * (gl1 + gl2) * j / (1.0 * min(gl1, gl2) * (1 + j))) ** (1.0 / k)


def dist_temp_func(cov, eps, k, l, cov_thres):
    if cov == "NA":
        return [1.0, 0]
    p = np.exp(-k * eps)
    copy_thres = int(1.0 * cov / cov_thres) + 1
    lam = 1.0 * cov * (l - k) / l
    if copy_thres == 1 or p == 1:
        return [1 - np.exp(-lam * p), lam * (1 - p)]
    else:
        s = [(lam * p) ** i / np.math.factorial(i) for i in range(copy_thres)]
        return [1 - np.exp(-lam * p) * sum(s), 0]


def estimate_dist(sample_1, sample_2, lib_1, lib_2, ce, le, ee, rl, k, cov_thres, tran):
    if sample_1 == sample_2 and lib_1 == lib_2:
        return sample_1, sample_2, 0.0
    sample_dir_1 = os.path.join(lib_1, sample_1)
    sample_dir_2 = os.path.join(lib_2, sample_2)
    msh_1 = os.path.join(sample_dir_1, sample_1 + ".msh")
    msh_2 = os.path.join(sample_dir_2, sample_2 + ".msh")
    dist_stderr = check_output(["mash", "dist", msh_1, msh_2], stderr=STDOUT, universal_newlines=True)
    j = float(dist_stderr.split()[4].split("/")[0]) / float(dist_stderr.split()[4].split("/")[1])
    gl_1 = le[sample_1]
    gl_2 = le[sample_2]
    if gl_1 == "NA" or gl_2 == "NA":
        gl_1 = 1
        gl_2 = 1
    cov_1 = ce[sample_1]
    cov_2 = ce[sample_2]
    eps_1 = ee[sample_1]
    eps_2 = ee[sample_2]
    l_1 = rl[sample_1]
    l_2 = rl[sample_2]
    r_1 = dist_temp_func(cov_1, eps_1, k, l_1, cov_thres)
    r_2 = dist_temp_func(cov_2, eps_2, k, l_2, cov_thres)
    wp = r_1[0] * r_2[0] * (gl_1 + gl_2) * 0.5
    zp = sum(r_1) * gl_1 + sum(r_2) * gl_2
    d = max(0, 1 - (1.0 * zp * j / (wp * (1 + j))) ** (1.0 / k))
    if tran:
        if d < 0.75:
            d = max(0, -0.75 * np.log(1 - 4.0 * d / 3.0))
        else:
            d = 5.0
    return sample_1, sample_2, d


def reference(args):
    # Creating a directory for reference library
    try:
        os.makedirs(args.l)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    # Creating a config file for references
    config_file = os.path.join(args.l, 'CONFIG')
    with open(config_file, mode='w') as f:
        f.write('kmer_length\t{0}\n'.format(args.k) + 'sketch_size\t{0}\n'.format(args.s) +
                'sketching_seed\t{0}\n'.format(args.S))

    # Making a list of sample names
    formats = ['.fq', '.fastq', '.fa', '.fna', '.fasta']
    files_names = [f for f in os.listdir(args.input_dir)
                   if True in (fnmatch.fnmatch(f, '*' + form) for form in formats)]
    samples_names = [f.rsplit('.f', 1)[0] for f in files_names]

    # Check if refs have duplicate entry
    if len(samples_names) != len(set(samples_names)):
        raise ValueError('Duplicate inputs (possibly same name with different extensions), please change '
                         'the file name(s) and try again')

    # Making a list of genome-skim files
    sequences = [os.path.join(args.input_dir, f) for f in files_names]

    # Initializing distance dataframe
    index = pd.MultiIndex.from_product([samples_names, samples_names], names=['sample', 'sample_2'])
    result_df = pd.DataFrame(columns=index)

    # Initializing coverage, genome length, error rate, and read length dictionaries
    cov_est = dict()
    len_est = dict()
    err_est = dict()
    read_len = dict()

    # Number of pools and threads for multi-processing
    n_pool = min(args.p, len(sequences))
    n_thread_cov = int(args.p / n_pool)
    n_proc_cov = n_pool * n_thread_cov
    n_pool_dist = min(args.p, len(sequences) ** 2)

    # Computing coverage, genome length, error rate, and read length
    sys.stderr.write('[skmer] Estimating coverages using {0} processors...\n'.format(n_proc_cov))
    pool_cov = mp.Pool(n_pool)
    results_cov = [pool_cov.apply_async(estimate_cov, args=(seq, args.l, args.k, args.e, n_thread_cov))
                   for seq in sequences]
    for result in results_cov:
        (name, coverage, genome_length, error_rate, read_length) = result.get(9999999)
        cov_est[name] = coverage
        len_est[name] = genome_length
        err_est[name] = error_rate
        read_len[name] = read_length
    pool_cov.close()
    pool_cov.join()
    
    # Sketching genome-skims
    sys.stderr.write('[skmer] Sketching sequences using {0} processors...\n'.format(n_pool))
    pool_sketch = mp.Pool(n_pool)
    results_sketch = [pool_sketch.apply_async(sketch, args=(seq, args.l, cov_est, err_est, args.k, args.s,
                                                            coverage_threshold, args.S)) for seq in sequences]
    for result in results_sketch:
        result.get(9999999)
    pool_sketch.close()
    pool_sketch.join()

    # Estimating pair-wise distances
    sys.stderr.write('[skmer] Estimating distances using {0} processors...\n'.format(n_pool_dist))
    pool_dist = mp.Pool(n_pool_dist)
    results_dist = [pool_dist.apply_async(estimate_dist, args=(s1, s2, args.l, args.l, cov_est, len_est, err_est,
                                                               read_len, args.k, coverage_threshold, args.t))
                    for s1 in samples_names for s2 in samples_names]

    for result in results_dist:
        dist_output = result.get(9999999)
        result_df[(dist_output[0], dist_output[1])] = [repr(dist_output[2])]

    # Writing distances to file
    sys.stderr.write('[skmer] Writing to file...\n')
    result_dfm = pd.melt(result_df, value_name='distance')
    result_mat = result_dfm.pivot(index='sample', columns='sample_2', values='distance')
    result_mat.to_csv(args.o + ".txt", sep='\t', mode='w')


def subsample(args):

    # Creating a directory for subsample
    try:
        os.makedirs(args.sub)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    # Making a list of sample names
    formats = ['.fq', '.fastq', '.fa', '.fna', '.fasta']
    files_names = [f for f in os.listdir(args.input_dir)
                   if True in (fnmatch.fnmatch(f, '*' + form) for form in formats)]
    samples_names = [f.rsplit('.f', 1)[0] for f in files_names]

    # Check if refs have duplicate entry
    if len(samples_names) != len(set(samples_names)):
        raise ValueError('Duplicate inputs (possibly same name with different extensions), please change '
                         'the file name(s) and try again')

    # Making a list of genome-skim files
    sequences = [os.path.join(args.input_dir, f) for f in files_names]

    # Initializing distance dataframe
    index = pd.MultiIndex.from_product([samples_names, samples_names], names=['sample', 'sample_2'])
    result_df = pd.DataFrame(columns=index)


    # Initializing coverage, genome length, error rate, and read length dictionaries
    cov_est = dict()
    len_est = dict()
    err_est = dict()
    read_len = dict()
    bs_kmer_sum = dict()
    sample_read_cnt = dict()

    # Number of pools and threads for multi-processing
    n_pool = min(args.p, len(sequences))
    n_thread_cov = int(args.p / n_pool)
    n_proc_cov = n_pool * n_thread_cov
    n_pool_dist = min(args.p, len(sequences) ** 2)

    # Computing coverage, genome length, error rate, read length and k-mer count
    sys.stderr.write('[skmer] Starting subsampling using {0} processors...\n'.format(n_proc_cov))
    pool_cov = mp.Pool(n_pool)
    results_cov = [pool_cov.apply_async(estimate_stats, args=(seq, n_thread_cov))
                   for seq in sequences]
    for result in results_cov:
        (name, coverage, genome_length, error_rate, read_length, rd_cnt) = result.get(9999999)
        cov_est[name] = coverage
        len_est[name] = genome_length
        err_est[name] = error_rate
        read_len[name] = read_length
        sample_read_cnt[name] = rd_cnt
        bs_kmer_sum[name] = args.s
    pool_cov.close()
    pool_cov.join()
    #print(sample_read_cnt)


    # Check whether inputs are reads or asseblies
    if "NA" in list(read_len.values()):
        input_data = 'assemblies'
    else:
        input_data = 'reads'

    
    ### Choose procedure for reads or assemblies ###
    sys.stderr.write('[skmer] Input processed as {}...\n'.format(input_data))

    # Compute block size

    np.random.seed(args.S)
    rand_seed_list = list(np.random.randint(low = 0, high = 4294967294, size = args.b))
    #print(rand_seed_list)

    bs_block_sz = {}
    bs_sample_sz = {}
    coef = args.c
    asm_sketch_sz = 0

    if input_data == 'reads':
       bs_sample_sz = sample_read_cnt
       for key, value in sample_read_cnt.items():
           bs_block_sz [key] = round((value)**(coef))
    else:
        bs_sample_sz = bs_kmer_sum
        mean_bs_kmer_count = np.mean(list(bs_kmer_sum.values()))
        asm_sketch_sz = round((mean_bs_kmer_count)**(coef))
        for key, value in bs_kmer_sum.items():
           bs_block_sz[key] = asm_sketch_sz

    #print(bs_block_sz)
    #print(bs_sample_sz)



    # Computing replicates
    for b in range (0, args.b):

        sys.stderr.write('[skmer] Computing replicate {0} using {1} processors...\n'.format(b, n_pool))

        # Creating replicate directory
        sub_rep = os.path.join(args.sub, "rep" + str(args.i))
        args.i +=1
        try:
            os.makedirs(sub_rep)
        except OSError as Error:
            if Error.errno != errno.EEXIST:
                raise

        # Creating replicate/library directory
        sub_lib = os.path.join(sub_rep, 'library')
        try:
            os.makedirs(sub_lib)
        except OSError as Error:
            if Error.errno != errno.EEXIST:
                raise

        # Update paths for subsampled replicate
        bs_sequences = [os.path.join(sub_rep, os.path.split(seq)[-1]) for seq in sequences]


        # Creating a config file for subsample  replicate
        config_file = os.path.join(sub_rep, 'CONFIG')
        with open(config_file, mode='w') as f:
            f.write('kmer_length\t{0}\n'.format(args.k) + 'sketch_size\t{0}\n'.format(args.s) +
                'sketching_seed\t{0}\n'.format(rand_seed_list[b]))


        # Write sample size dictionary to file
        np.save(os.path.join(sub_rep, 'block_size.npy'), bs_block_sz)        
        np.save(os.path.join(sub_rep, 'sample_size.npy'), bs_sample_sz)

        # Update  coverage and error estimates for subsample
        if input_data == 'reads':


            # Generate subsample replicates and save to bootstrap directory
            pool_sketch = mp.Pool(n_pool)
            results_sketch = [pool_sketch.apply_async(sample_reads, args=(seq, rand_seed_list[b], bs_block_sz[(os.path.split(seq)[-1]).rsplit('.f', 1)[0]], sub_rep)) for seq in sequences]
            for result in results_sketch:
                result.get(9999999)
            pool_sketch.close()
            pool_sketch.join()


            # Computing coverage, genome length, error rate, and read length of replicates  using reference function
            pool_cov = mp.Pool(n_pool)
            results_cov = [pool_cov.apply_async(estimate_cov, args=(seq, sub_lib, args.k, args.e, n_thread_cov))
                       for seq in bs_sequences]
            for result in results_cov:
                (name, coverage, genome_length, error_rate, read_length) = result.get(9999999)
                cov_est[name] = coverage
                len_est[name] = genome_length
                err_est[name] = error_rate
                read_len[name] = read_length
            pool_cov.close()
            pool_cov.join()


            # Sketching genome-skims
            pool_sketch = mp.Pool(n_pool)
            #reads_sketch_sz = 100000
            results_sketch = [pool_sketch.apply_async(sketch, args=(seq, sub_lib, cov_est, err_est, args.k, args.s,
                                                            coverage_threshold, rand_seed_list[b])) for seq in bs_sequences]
            for result in results_sketch:
                result.get(9999999)
            pool_sketch.close()
            pool_sketch.join()


            # Estimating pair-wise distances
            pool_dist = mp.Pool(n_pool_dist)
            results_dist = [pool_dist.apply_async(estimate_dist, args=(s1, s2, sub_lib, sub_lib, cov_est, len_est, err_est,
                                                                   read_len, args.k, coverage_threshold, args.t))
                        for s1 in samples_names for s2 in samples_names]

            for result in results_dist:
                dist_output = result.get(9999999)
                result_df[(dist_output[0], dist_output[1])] = [repr(dist_output[2])]


        else:

            # Prepare genome-skims directory structure
            pool_sketch = mp.Pool(n_pool)
            results_sketch = [pool_sketch.apply_async(create_sketch_dir, args=(seq, sub_lib, cov_est, len_est, err_est, 
                                                                               read_len, args.t)) for seq in sequences]
            pool_sketch.close()
            pool_sketch.join()



            # Sketching genome-skims
            pool_sketch = mp.Pool(n_pool)
            results_sketch = [pool_sketch.apply_async(sketch, args=(seq, sub_lib, cov_est, err_est, args.k, asm_sketch_sz,
                                                            coverage_threshold, rand_seed_list[b])) for seq in sequences]
            for result in results_sketch:
                result.get(9999999)
            pool_sketch.close()
            pool_sketch.join()



            # Estimating pair-wise distances
            pool_dist = mp.Pool(n_pool_dist)
            results_dist = [pool_dist.apply_async(estimate_dist, args=(s1, s2, sub_lib, sub_lib, cov_est, len_est, err_est,
                                                                   read_len, args.k, coverage_threshold, args.t))
                        for s1 in samples_names for s2 in samples_names]

            for result in results_dist:
                dist_output = result.get(9999999)
                result_df[(dist_output[0], dist_output[1])] = [repr(dist_output[2])]



        # Writing distances to file
        sys.stderr.write('[skmer] Writing to file...\n')
        result_dfm = pd.melt(result_df, value_name='distance')
        result_mat = result_dfm.pivot(index='sample', columns='sample_2', values='distance')
        final_path = os.path.join(sub_rep, "dimtrx_rep" + ".txt")
        result_mat.to_csv(final_path, float_format='%f', sep='\t', mode='w')


        # Cleaning up 
        if args.fa:
                for fi in bs_sequences:
                    try:
                        os.remove(fi)
                    except OSError:
                        pass

        if args.msh:
            sketch_fi = []
            for (dirpath, dirnames, filenames) in os.walk(sub_lib):
                sketch_fi += [os.path.join(dirpath, file) for file in filenames if file.endswith(".msh") ]
            #print(sketch_fi)
            for fi in sketch_fi:
                try:
                    os.remove(fi)
                except OSError:
                    pass


    # Clean up subsample folders 
    #import shutil
    #shutil.rmtree(sub_lib)
    #shutil.rmtree(args.bs)
    #shutil.rmtree(args.l)


def correction(args):
     


    # Making a list of sample names
    try:
        #with open(args.main,"r") as f:
        df = pd.read_csv(args.main, header = 0, sep='\t', skiprows = 0)
        samples_names = list(df.iloc[:,0])
        #print(samples_names)
    except:
        raise ValueError('Please check file name for main distance matrix and try again')


    # Initializing distance dataframe
    index = pd.MultiIndex.from_product([samples_names, samples_names], names=['sample', 'sample_2'])
    result_df = pd.DataFrame(columns=index)
    combo_result_df = pd.DataFrame()


    # Round distances up to 12 digits since fastme doesn't except more than 12 decimals
    no_strap_dfm = pd.melt(df, id_vars=['sample'], value_vars= list(df.columns[1:]) )
    no_strap_dfm.rename(columns={'variable':'sample_2'}, inplace=True) 
    no_strap_dfm.rename(columns={'value':'no_strapped_dist'}, inplace=True)    

    decimals = 12
    no_strap_dfm['no_strapped_dist'] = no_strap_dfm['no_strapped_dist'].apply(pd.to_numeric, errors='coerce')
    no_strap_dfm['no_strapped_dist'] =  no_strap_dfm['no_strapped_dist'].apply(lambda x: round(x, decimals))

    no_strap_mat = no_strap_dfm.pivot(index='sample', columns='sample_2', values='no_strapped_dist')
    no_strap_mat.to_csv(os.path.splitext(args.main)[0] + "_cor_" +  ".txt", float_format='%f', sep='\t', mode='w')

    # List replicate directories
    try:
        for dir in [name for name in os.listdir(args.sub) if 'rep' in name]:
            print(dir)
            rep_mtrx = os.path.join(args.sub, dir, "dimtrx_rep.txt")
            df = pd.read_csv(rep_mtrx, header = 0, sep='\t', skiprows = 0)

            # Append estimates to combo dataframe
            result_dfm = pd.melt(df, id_vars=['sample'], value_vars= list(df.columns[1:]) )
            result_dfm.rename(columns={'variable':'sample_2'}, inplace=True)
            result_dfm.rename(columns={'value':'uncorrected_dist'}, inplace=True)
            result_dfm['rep'] = int(dir.split('rep', 1)[-1])

            # Load dictionaries
            bs_block_sz = np.load(os.path.join(args.sub, dir, 'block_size.npy'), allow_pickle='TRUE').item()
            bs_sample_sz = np.load(os.path.join(args.sub, dir, 'sample_size.npy'), allow_pickle='TRUE').item()
            
            result_dfm['b_s1'] = result_dfm['sample'].map(bs_block_sz)
            result_dfm['b_s2'] = result_dfm['sample_2'].map(bs_block_sz)
            result_dfm['b_s1'] = result_dfm['b_s1'].apply(pd.to_numeric, errors='coerce')
            result_dfm['b_s2'] = result_dfm['b_s2'].apply(pd.to_numeric, errors='coerce')
            result_dfm['b_mean'] = result_dfm[['b_s1', 'b_s2']].mean(axis=1, skipna=True)
            
            result_dfm['n_s1'] = result_dfm['sample'].map(bs_sample_sz)
            result_dfm['n_s2'] = result_dfm['sample_2'].map(bs_sample_sz)
            result_dfm['n_s1'] = result_dfm['n_s1'].apply(pd.to_numeric, errors='coerce')
            result_dfm['n_s2'] = result_dfm['n_s2'].apply(pd.to_numeric, errors='coerce')
            result_dfm['N_mean'] = result_dfm[['n_s1', 'n_s2']].mean(axis=1, skipna=True)
            combo_result_df = combo_result_df.append(result_dfm, ignore_index = True)
            #print(combo_result_df)

    except:
        raise ValueError('Please check subsample directory and try again')


   
    # Computing distance correction

    combo_result_df['uncorrected_dist'] = combo_result_df['uncorrected_dist'].apply(pd.to_numeric, errors='coerce')
    res = combo_result_df.groupby(['sample', 'sample_2'], as_index=False)['uncorrected_dist'].mean()
    res.rename({'uncorrected_dist': 'subsample_mean_dist'}, axis=1, inplace=True)    
    new_df = pd.merge(combo_result_df, res,  how='left', left_on=['sample','sample_2'], right_on = ['sample','sample_2'])

    new_df_out = pd.merge(new_df, no_strap_dfm,  how='left', left_on=['sample','sample_2'], right_on = ['sample','sample_2'])
    new_df_out['no_strapped_dist'] = new_df_out['no_strapped_dist'].apply(pd.to_numeric, errors='coerce')
    
    new_df_out['corrected_dist'] = ((new_df_out['b_mean']/new_df_out['N_mean'])**(1/2))*(new_df_out['uncorrected_dist']-new_df_out['subsample_mean_dist'])+new_df_out['no_strapped_dist']
    new_df_out['corrected_dist_cons'] = ((new_df_out['b_mean']/new_df_out['N_mean'])**(1/2))*(new_df_out['uncorrected_dist']-new_df_out['subsample_mean_dist'])+new_df_out['subsample_mean_dist']
    
    #replace negative values with 0.0 so fastme can handle matrices
    new_df_out.corrected_dist = np.where(new_df_out.corrected_dist < 0, 0.0, new_df_out.corrected_dist)
    new_df_out.corrected_dist_cons = np.where(new_df_out.corrected_dist_cons < 0, 0.0, new_df_out.corrected_dist_cons)

    # round distances up to 12 digits since fastme doesn't except more than 12 decimals
    new_df_out['corrected_dist'] =  new_df_out['corrected_dist'].apply(lambda x: round(x, decimals))
    new_df_out['corrected_dist_cons'] =  new_df_out['corrected_dist_cons'].apply(lambda x: round(x, decimals))
    new_df_out.to_csv(os.path.join(args.sub, "_summary" + ".csv"), float_format='%f', sep=',', mode='w')



    # Writing distances to file
    sys.stderr.write('[skmer] Writing to file...\n')
    
    b_list = list((new_df_out.loc[:, 'rep']).unique())
    #print(b_list)
    for b in b_list:
        sub_dfm = new_df_out.loc[(new_df_out['rep'] == b)]
        
        #-mean+main
        sub_dfm_main = sub_dfm[['sample','sample_2', 'corrected_dist']]
        result_mat_main = sub_dfm_main.pivot(index='sample', columns='sample_2', values='corrected_dist')
        result_mat_main.to_csv(os.path.join(args.sub, "rep" + str(b),  "dimtrx_rep_cor.txt"), float_format='%f', sep='\t', mode='w')
        
        #-mean+mean
        sub_dfm_cons = sub_dfm[['sample','sample_2', 'corrected_dist_cons']]
        result_mat_cons = sub_dfm_cons.pivot(index='sample', columns='sample_2', values='corrected_dist_cons')
        result_mat_cons.to_csv(os.path.join(args.sub, "rep" + str(b),  "dimtrx_rep_cor_cons.txt"), float_format='%f', sep='\t', mode='w')


def distance(args):
    # Loading reference config
    config_file = os.path.join(args.library, 'CONFIG')
    with open(config_file) as f:
        config = f.read()
    kl = int(config.split('\n')[0].split('\t')[1])

    # Making a list of reference samples
    refs = [item for item in os.listdir(args.library) if os.path.isdir(os.path.join(args.library, item))]

    # Initializing distance dataframe
    index = pd.MultiIndex.from_product([refs, refs], names=['sample', 'sample_2'])
    result_df = pd.DataFrame(columns=index)

    # Loading coverage, genome length, error rate, and read length information
    cov_est = dict()
    len_est = dict()
    err_est = dict()
    read_len = dict()
    for ref in refs:
        ref_dir = os.path.join(args.library, ref)
        info_file = os.path.join(ref_dir, ref + '.dat')
        with open(info_file) as f:
            info = f.read()
        cov_value = info.split('\n')[0].split('\t')[1]
        gl_value = info.split('\n')[1].split('\t')[1]
        if cov_value == "NA":
            if gl_value == "NA":
                cov_est[ref] = "NA"
                len_est[ref] = "NA"
                err_est[ref] = "NA"
                read_len[ref] = int(info.split('\n')[3].split('\t')[1])
            else:
                cov_est[ref] = "NA"
                len_est[ref] = int(info.split('\n')[1].split('\t')[1])
                err_est[ref] = 0
                read_len[ref] = "NA"
        else:
            cov_est[ref] = float(info.split('\n')[0].split('\t')[1])
            len_est[ref] = int(info.split('\n')[1].split('\t')[1])
            err_est[ref] = float(info.split('\n')[2].split('\t')[1])
            read_len[ref] = int(info.split('\n')[3].split('\t')[1])

    # Number of pools and threads for multi-processing
    n_pool_dist = min(args.p, len(refs) ** 2)

    # Estimating pair-wise distances
    sys.stderr.write('[skmer] Estimating distances using {0} processors...\n'.format(n_pool_dist))
    pool_dist = mp.Pool(n_pool_dist)
    results_dist = [pool_dist.apply_async(estimate_dist, args=(r1, r2, args.library, args.library, cov_est, len_est,
                                                               err_est, read_len, kl, coverage_threshold, args.t))
                    for r1 in refs for r2 in refs]

    for result in results_dist:
        dist_output = result.get(9999999)
        result_df[(dist_output[0], dist_output[1])] = [repr(dist_output[2])]

    # Writing distances to file
    sys.stderr.write('[skmer] Writing to file...\n')
    result_dfm = pd.melt(result_df, value_name='distance')
    result_mat = result_dfm.pivot(index='sample', columns='sample_2', values='distance')
    result_mat.to_csv(args.o + ".txt", sep='\t', mode='w')


def query(args):
    # Loading reference config
    config_file = os.path.join(args.library, 'CONFIG')
    with open(config_file) as f:
        config = f.read()
    kl = int(config.split('\n')[0].split('\t')[1])
    ss = int(config.split('\n')[1].split('\t')[1])
    seed = int(config.split('\n')[2].split('\t')[1])

    # Creating a directory for the query
    sample = os.path.basename(args.input).rsplit('.f', 1)[0]
    sample_dir = os.path.join(os.getcwd(), sample)
    try:
        os.makedirs(sample_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    # Making a list of references samples
    refs = [item for item in os.listdir(args.library) if os.path.isdir(os.path.join(args.library, item))]

    # Check if the sample is already in the refs
    if sample in set(refs):
        raise ValueError('A reference sample exists with the same name as the query, please change '
                         'the name of query file {0} and try again'.format(sample))

    # Initializing distances series
    result_s = pd.Series(index=refs, name=sample)

    # Loading coverage, genome length, error rate, and read length information
    cov_est = dict()
    len_est = dict()
    err_est = dict()
    read_len = dict()
    for ref in refs:
        ref_dir = os.path.join(args.library, ref)
        info_file = os.path.join(ref_dir, ref + '.dat')
        with open(info_file) as f:
            info = f.read()
        cov_value = info.split('\n')[0].split('\t')[1]
        gl_value = info.split('\n')[1].split('\t')[1]
        if cov_value == "NA":
            if gl_value == "NA":
                cov_est[ref] = "NA"
                len_est[ref] = "NA"
                err_est[ref] = "NA"
                read_len[ref] = int(info.split('\n')[3].split('\t')[1])
            else:
                cov_est[ref] = "NA"
                len_est[ref] = int(info.split('\n')[1].split('\t')[1])
                err_est[ref] = 0
                read_len[ref] = "NA"
        else:
            cov_est[ref] = float(info.split('\n')[0].split('\t')[1])
            len_est[ref] = int(info.split('\n')[1].split('\t')[1])
            err_est[ref] = float(info.split('\n')[2].split('\t')[1])
            read_len[ref] = int(info.split('\n')[3].split('\t')[1])

    # Number of pools for multi-processing
    n_pool_dist = min(args.p, len(refs))

    # Computing the coverage, genome length, error rate, and read length of query sample
    sys.stderr.write('[skmer] Estimating the coverage using {0} processors...\n'.format(args.p))
    (dummy, coverage, genome_length, error_rate, read_length) = estimate_cov(args.input, os.getcwd(), kl, args.e,
                                                                             args.p)
    cov_est[sample] = coverage
    len_est[sample] = genome_length
    err_est[sample] = error_rate
    read_len[sample] = read_length

    # Sketching the query genome-skim
    sys.stderr.write('[skmer] Sketching the genome-skim...\n')
    sketch(args.input, os.getcwd(), cov_est, err_est, kl, ss, coverage_threshold, seed)

    # Estimating pair-wise distances
    sys.stderr.write('[skmer] Estimating distances using {0} processors...\n'.format(n_pool_dist))
    pool_dist = mp.Pool(n_pool_dist)
    results_dist = [pool_dist.apply_async(estimate_dist, args=(sample, ref, os.getcwd(), args.library, cov_est, len_est,
                                                               err_est, read_len, kl, coverage_threshold, args.t))
                    for ref in refs]
    for result in results_dist:
        dist_output = result.get(9999999)
        result_s[dist_output[1]] = dist_output[2]

    # Writing distances to file
    sys.stderr.write('[skmer] Writing to file...\n')
    result_s.sort_values(inplace=True)
    result_sr = result_s.apply(repr)
    result_sr.to_csv('{0}-{1}.txt'.format(args.o, sample.lower()), sep='\t', mode='w')

    # Adding query to the reference library
    if args.a:
        try:
            shutil.copytree(sample_dir, os.path.join(args.library, sample))
        except shutil.Error as e:
            print('Directory not copied. Error: %s' % e)
        except OSError as e:
            print('Directory not copied. Error: %s' % e)

    shutil.rmtree(sample_dir)


def main():
    # Input arguments parser
    parser = argparse.ArgumentParser(description='{0} - Estimating genomic distances between '.format(__version__) +
                                                 'genome-skims',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument('-v', '--version', action='store_true', help='print the current version')
    parser.add_argument('--debug', action='store_true', help='Print the traceback when an exception is raised')
    subparsers = parser.add_subparsers(title='commands',
                                       description='reference   Process a library of reference genome-skims or assemblies\n'
                                                   'distance    Compute pairwise distances for a processed library\n'
                                                   'query       Compare a genome-skim or assembly against a reference library\n'
                                                   'subsample   Performs  subsample on a library of reference genome-skims or assemblies\n'
                                                   'correct     Performs correction of subsampled distance matrices obtained for reference' 
                                                   ' genome-skims or assemblies'
                                                   ,
                                       help='Run skmer {commands} [-h] for additional help',
                                       dest='{commands}')

    # To make sure that subcommand is required in python >= 3.3
    python_version = sys.version_info
    if (python_version[0] * 10 + python_version[1]) >= 33:
        subparsers.required = True

    # Reference command subparser
    parser_ref = subparsers.add_parser('reference',
                                       description='Process a library of reference genome-skims or assemblies')
    parser_ref.add_argument('input_dir',
                            help='Directory of input genome-skims or assemblies (dir of .fastq/.fq/.fa/.fna/.fasta files)')
    parser_ref.add_argument('-l', default=os.path.join(os.getcwd(), 'library'),
                            help='Directory of output (reference) library. Default: working_directory/library')
    parser_ref.add_argument('-o', default='ref-dist-mat',
                            help='Output (distances) prefix. Default: ref-dist-mat')
    parser_ref.add_argument('-k', type=int, choices=list(range(1, 32)), default=31, help='K-mer length [1-31]. ' +
                                                                                         'Default: 31', metavar='K')
    parser_ref.add_argument('-s', type=int, default=10 ** 5, help='Sketch size. Default: 100000')
    parser_ref.add_argument('-S', type=int, default=42, help='Sketching random seed. Default: 42')
    parser_ref.add_argument('-e', type=float, help='Base error rate. By default, the error rate is automatically '
                                                   'estimated.')
    parser_ref.add_argument('-t', action='store_true',
                            help='Apply Jukes-Cantor transformation to distances. Output 5.0 if not applicable')
    parser_ref.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_ref.set_defaults(func=reference)

    # Subsample command subparser
    parser_bt = subparsers.add_parser('subsample',
                                       description='Performs  subsample on a library of reference genome-skims or assemblies')
    parser_bt.add_argument('input_dir',
                            help='Directory of input genome-skims or assemblies (dir of .fastq/.fq/.fa/.fna/.fasta files)')
    #parser_bt.add_argument('-l', default=os.path.join(os.getcwd(), 'library'),
    #                        help='Directory of output (reference) library. Default: working_directory/library')
    parser_bt.add_argument('-sub', default=os.path.join(os.getcwd(), 'subsample'),
                            help='Directory of output for subsample replicates. Default: working_directory/subsample')
    parser_bt.add_argument('-fa', action='store_false',
                            help='Save subsampled genome-skims. Default: false')
    parser_bt.add_argument('-msh', action='store_false',
                            help='Save sketches. Default: false')
   # parser_bt.add_argument('-o', default='ref-dist-mat',
   #                         help='Output (distances) prefix. Default: ref-dist-mat')
    parser_bt.add_argument('-k', type=int, choices=list(range(1, 32)), default=31, help='K-mer length [1-31]. ' +
                                                                                         'Default: 31', metavar='K')
    parser_bt.add_argument('-s', type=int, default=10 ** 5, help='Sketch size. Default: 100000')
    parser_bt.add_argument('-S', type=int, default=42, help='Sketching random seed. Default: 42')
    parser_bt.add_argument('-i', type=int, default=0, help='Start index of subsampled replicate (eg 5 for dir rep5). Default: 0')
    parser_bt.add_argument('-b', type=int, default=100, help='Number of subsampled replicates. Default: 100')    
    parser_bt.add_argument('-c', type=float, default=0.9, help='Exponent value for subsampling. Default: 0.9')
    parser_bt.add_argument('-e', type=float, help='Base error rate. By default, the error rate is automatically '
                                                   'estimated.')
    parser_bt.add_argument('-t', action='store_true',
                            help='Apply Jukes-Cantor transformation to distances. Output 5.0 if not applicable')
    parser_bt.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_bt.set_defaults(func=subsample)    


    # Correction command subparser
    parser_cor = subparsers.add_parser('correct',
                                       description='Performs correction of subsampled distance matrices obtained for reference genome-skims or assemblies')
    parser_cor.add_argument('-main',
                            help='Distance matrix of main estimate')
    parser_cor.add_argument('-sub', default=os.path.join(os.getcwd(), 'subsample'),
                            help='Directory of output for subsample replicates. Default: working_directory/subsample')
    parser_cor.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_cor.set_defaults(func=correction)

    # Distance command subparser
    parser_dist = subparsers.add_parser('distance', description='Compute the distance matrix for a processed library')
    parser_dist.add_argument('library', help='Directory of the input (processed) library')
    parser_dist.add_argument('-o', default='ref-dist-mat',
                             help='Output (distances) prefix. Default: ref-dist-mat')
    parser_dist.add_argument('-t', action='store_true',
                             help='Apply Jukes-Cantor transformation to distances. Output 5.0 if not applicable')
    parser_dist.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                             help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                  'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_dist.set_defaults(func=distance)

    # query command subparser
    parser_qry = subparsers.add_parser('query',
                                       description='Compare an input genome-skim or assembly against a reference library')
    parser_qry.add_argument('input', help='Input (query) genome-skim or assembly (a .fastq/.fq/.fa/.fna/.fasta file)')
    parser_qry.add_argument('library', help='Directory of (reference) library')
    parser_qry.add_argument('-a', action='store_true',
                            help='Add the processed input (query) to the (reference) library')
    parser_qry.add_argument('-o', default='dist',
                            help='Output (distances) prefix. Default: dist')
    parser_qry.add_argument('-e', type=float, help='Base error rate. By default, the error rate is automatically '
                                                   'estimated.')
    parser_qry.add_argument('-t', action='store_true',
                            help='Apply Jukes-Cantor transformation to distances. Output 5.0 if not applicable')
    parser_qry.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_qry.set_defaults(func=query)

    args = parser.parse_args()

    # Handling traceback on exceptions
    def exception_handler(exception_type, exception, traceback, debug_hook=sys.excepthook):
        if args.debug:
            debug_hook(exception_type, exception, traceback)
        else:
            print("{0}: {1}".format(exception_type.__name__, exception))

    sys.excepthook = exception_handler

    args.func(args)


if __name__ == "__main__":
    main()
