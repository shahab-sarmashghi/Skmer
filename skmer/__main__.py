#! /usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
from scipy.optimize import newton, curve_fit
import argparse
import os
import shutil
import fnmatch
import sys
import errno
import pandas as pd
from subprocess import call, check_output, STDOUT
import multiprocessing as mp

__version__ = 'skmer-pop v3'

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
    return int(round(1.0 * total_length / n_reads)), max_length, total_length


def cov_temp_func(x, r, p, k, l):
    lam = x * (1.0 * (l - k)) / l
    return lam * (p ** 2) * np.exp(-lam * p) - 2 * r * (p * np.exp(-lam * p) + 1 - p)


def poisson(mu, i):
    return np.exp(-mu) * mu ** i / np.math.factorial(i)


def hist_func_new(x, eps, q_1, q_2, q_3, q_4, q_5):
    q = np.array([q_1, q_2, q_3, q_4, q_5])
    J = np.linspace(1, len(q), len(q))
    M = np.zeros(len(x))
    lam = x[0][1] / np.dot(q, J)
    zeta = lam * (1 - eps) ** x[0][0]
    for i in range(1, len(x)):
        if i==1:
            P = [poisson(j * zeta, 1) + (lam - zeta) * j for j in J]
            M[i] = np.dot(q, P)
        else:
            P = [poisson(j * zeta, x[i]) for j in J]
            M[i] = np.dot(q, P)
    return M[1:]


def estimate_cov(sequence, lib, k, e, nth):
    sample = os.path.basename(sequence).rsplit('.f', 1)[0]
    sample_dir = os.path.join(lib, sample)
    try:
        os.makedirs(sample_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise
    info_file = os.path.join(sample_dir, sample + '.dat')

    (l, ml, tl) = sequence_stat(sequence)
    if ml > seq_len_threshold:
        cov = "NA"
        g_len = tl
        eps = 0
        l = "NA"
 ## Add repeats
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
    # with open(histo_file) as f:
    #     histo_stderr = f.read()
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
## Add repeats
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
        lam_0 = (np.exp(-gam) * (gam ** ind) / np.math.factorial(ind)) * count[1] / count[ind] + gam * (1 - np.exp(-gam))
        eps_0 = 1 - (gam / lam_0) ** (1.0 / k)
        cov_0 = (1.0 * l / (l - k)) * lam_0
        zeta_0 = lam_0 * (1 - eps_0) ** k
## Fix this for count and multiplicit
    tot_seq = 1.0 * ksum * l / (l - k)
    g_len_0 = int(tot_seq / cov_0)
    n_rep_term = 5
    n_bins = 20
    x_0 = list(g_len_0 * np.array([0.9, 0.09, 0.009, 0.0005, 0.005]))
    xdata_0 = list(range(1, n_bins + 1))
    ydata_0 = [count[x] for x in xdata_0]
    xdata_0 = [(k, 1.0 * ksum)] + xdata_0
    popt, pcov = curve_fit(hist_func_new, xdata_0, ydata_0, p0=[eps_0]+x_0,
                           bounds=(0, [10*eps_0, 10*g_len_0, 10*g_len_0, 10*g_len_0, 10*g_len_0, 10*g_len_0]))
    eps = popt[0]
    q = popt[1:]
    J = np.linspace(1, len(q), len(q))
    g_len = np.dot(q, J)
    lam = xdata_0[0][1] / g_len
    cov = (1.0 * l / (l - k)) * lam
    rep_prob = list(1.0 * q / sum(q))

    if eps > error_rate_threshold or eps < 0:
        cov = "NA"
        g_len = "NA"
        eps = "NA"
        with open(info_file, mode='w') as f:
            f.write('coverage\t{0}\n'.format(cov) + 'genome_length\t{0}\n'.format(g_len) +
                    'error_rate\t{0}\n'.format(eps) + 'read_length\t{0}\n'.format(l))
        return sample, cov, g_len, eps, l

    q_string = ''
    for prob in rep_prob:
        q_string += repr(prob) + '\t'
    with open(info_file, mode='w') as f:
        f.write('coverage\t{0}\n'.format(repr(cov)) + 'genome_length\t{0}\n'.format(g_len) +
                'error_rate\t{0}\n'.format(repr(eps)) + 'read_length\t{0}\n'.format(l) +
                'repeat_profile\t{0}\n'.format(q_string.strip()))

## Output repeats
    # return sample, cov, g_len, eps, l, rep_prob
    return sample, cov, g_len, eps, l


def poisson_percent(mu, p):
    P = 0
    for i in range(int(mu)+1):
        P += poisson(mu, i)
        if P > 1-p:
            return max(1, i)
    return 1


def sketch(sequence, lib, ce, ee, k, s, cov_thres):
    sample = os.path.basename(sequence).rsplit('.f', 1)[0]
    sample_dir = os.path.join(lib, sample)
    msh = os.path.join(sample_dir, sample)
    cov = ce[sample]
    eps = ee[sample]
    if cov == "NA" and eps == 0:
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-o", msh, sequence], stderr=open(os.devnull, 'w'))
        return
    elif eps == "NA":
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-r", "-o", msh, sequence], stderr=open(os.devnull, 'w'))
        return
    copy_thres = int(cov / cov_thres) + 1
    if cov < cov_thres or eps == 0.0:
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-r", "-o", msh, sequence], stderr=open(os.devnull, 'w'))
    else:
        call(["mash", "sketch", "-m", str(copy_thres), "-k", str(k), "-s", str(s), "-o", msh, sequence],
             stderr=open(os.devnull, 'w'))
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
        f.write('kmer_length\t{0}\n'.format(args.k) + 'sketch_size\t{0}\n'.format(args.s))

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
                                                            coverage_threshold)) for seq in sequences]
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
            len_est[ref] = int(float(info.split('\n')[1].split('\t')[1]))
            err_est[ref] = float(info.split('\n')[2].split('\t')[1])
            read_len[ref] = int(float(info.split('\n')[3].split('\t')[1]))

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
            len_est[ref] = int(float(info.split('\n')[1].split('\t')[1]))
            err_est[ref] = float(info.split('\n')[2].split('\t')[1])
            read_len[ref] = int(float(info.split('\n')[3].split('\t')[1]))

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
    sketch(args.input, os.getcwd(), cov_est, err_est, kl, ss, coverage_threshold)

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
        os.rename(sample_dir, os.path.join(args.library, sample))
    # To save processed query by default
    else:
        os.rename(sample_dir, os.path.join(args.O, sample))


def main():
    # Input arguments parser
    parser = argparse.ArgumentParser(description='{0} - Estimating gonomic distances between '.format(__version__) +
                                                 'genome-skims',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument('-v', '--version', action='store_true', help='print the current version')
    parser.add_argument('--debug', action='store_true', help='Print the traceback when an exception is raised')
    subparsers = parser.add_subparsers(title='commands',
                                       description='reference   Process a library of reference genome-skims or assemblies\n'
                                                   'distance    Compute pairwise distances for a processed library\n'
                                                   'query       Compare a genome-skim or assembly against a reference '
                                                   'library',
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
    parser_ref.add_argument('-s', type=int, default=10 ** 7, help='Sketch size. Default: 10000000')
    parser_ref.add_argument('-e', type=float, help='Base error rate. By default, the error rate is automatically '
                                                   'estimated.')
    parser_ref.add_argument('-t', action='store_true',
                            help='Apply Jukes-Cantor transformation to distances. Output 5.0 if not applicable')
    parser_ref.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_ref.set_defaults(func=reference)

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
    parser_qry.add_argument('-O', default=os.getcwd(),
                            help='Output (processed) query path. Default: working_directory/')
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
