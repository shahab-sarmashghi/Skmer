#! /usr/bin/env python
# -*- coding: utf-8 -*-


# __version__ = 'skmer2'


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


def compute_read_length(skim):
    with open(skim) as f:
        for line in f:
            if line.startswith('@'):
                return len(next(f).strip())
            elif line.startswith('>'):
                l = len(next(f).strip())
                next_line = next(f)
                while not next_line.startswith('>'):
                    l += len(next_line.strip())
                    next_line = next(f)
                return l


def cov_func(x, r, p, k, l):
    lam = x * (1.0 * (l-k)) / l
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


def estimate_cov(skim, lib, k, e, nth):
    sample = os.path.basename(skim).split('.f')[0]
    sample_dir = os.path.join(lib, sample)
    try:
        os.makedirs(sample_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise
    mercnt = os.path.join(sample_dir, sample + '.jf')
    histo_file = os.path.join(sample_dir, sample + '.hist')
    call(["jellyfish", "count", "-m", str(k), "-s", "100M", "-t", str(nth), "-C", "-o", mercnt, skim],
         stderr=open(os.devnull, 'w'))
    histo_stderr = check_output(["jellyfish", "histo", "-h", "1000000", mercnt], stderr=STDOUT, universal_newlines=True)
    with open(histo_file, mode='w') as f:
        f.write(histo_stderr)
    os.remove(mercnt)
    count = {0: 0}
    multiplicity = [0]
    ksum = 0
    for item in histo_stderr.split('\n')[:-1]:
        count[int(item.split()[0])] = int(item.split()[1])
        multiplicity.append(float(item.split()[0]))
        ksum += int(item.split()[1]) * int(item.split()[0])
    if len(count) < 3:
        raise ValueError('Coverage of {0} is too low; unable to estimate it'.format(sample))
    count_ef = {x: count[x] for x in count if x not in [0, 1]}
    ind = min(max(count_ef, key=count_ef.get), len(count) - 2)  # The second argument needs to be safeguarded

    l = compute_read_length(skim)
    if e is not None:
        pass
        # Needs to be completed later
        # eps = e
        # p0 = np.exp(-k * eps)
        # if ind < 2:
        #     r21 = 1.0 * count[2] / count[1]
        #     cov = newton(cov_func, 0.05, args=(r21, p0, k, l))
        # else:
        #     cov = (1.0 / p0) * (1.0 * l / (l - k)) * (ind + 1) * count[ind + 1] / count[ind]
    elif ind < 2:
        raise ValueError('Not enough information to co-estimate coverage and error rate of {0}'.format(sample))
    else:
        gam = 1.0 * (ind + 1) * count[ind + 1] / count[ind]
        lam_0 = (np.exp(-gam) * (gam ** ind) / np.math.factorial(ind)) * count[1] / count[ind] + gam * (1 - np.exp(-gam))
        eps_0 = 1 - (gam / lam_0) ** (1.0 / k)
        cov_0 = (1.0 * l / (l - k)) * lam_0
        zeta_0 = lam_0 * (1 - eps_0) ** k

    tot_seq = 1.0 * ksum * l / (l - k)
    g_len_0 = int(tot_seq / cov_0)
    n_rep_term = 5
    n_bins = 20
    x_0 = list(g_len_0 * np.array([0.9, 0.09, 0.009, 0.0005, 0.005]))
    xdata_0 = multiplicity[1:n_bins + 1]
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

    info_file = os.path.join(sample_dir, sample + '.dat')
    q_string = ''
    for prob in rep_prob:
        q_string += repr(prob) + '\t'
    with open(info_file, mode='w') as f:
        f.write('coverage\t{0}\n'.format(repr(cov)) + 'genome_length\t{0}\n'.format(g_len) +
                'error_rate\t{0}\n'.format(repr(eps)) + 'read_length\t{0}\n'.format(l) +
                'repeat_profile\t{0}\n'.format(q_string.strip()))

    # with open(info_file) as f:
    #     info = f.read()
    # cov = float(info.split('\n')[0].split('\t')[1])
    # g_len = float(info.split('\n')[1].split('\t')[1])
    # eps = float(info.split('\n')[2].split('\t')[1])
    # l = float(info.split('\n')[3].split('\t')[1])
    # rep_prob = [float(x) for x in info.split('\n')[4].split('\t')[1:]]

    return sample, cov, g_len, eps, l, rep_prob


def poisson_percent(mu, p):
    P = 0
    for i in range(int(mu)+1):
        P += poisson(mu, i)
        if P > 1-p:
            return max(1, i)
    return 1


def sketch(skim, lib, ce, ee, k, s, rl):
    sample = os.path.basename(skim).split('.f')[0]
    sample_dir = os.path.join(lib, sample)
    msh = os.path.join(sample_dir, sample)
    cov = ce[sample]
    eps = ee[sample]
    l = rl[sample]
    p = np.exp(-k * eps)
    lam = cov * (l - k) / l
    copy_thres = poisson_percent(lam * p, 0.95)
    if eps == 0.0:
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-r", "-o", msh, skim], stderr=open(os.devnull, 'w'))
    else:
        call(["mash", "sketch", "-m", str(copy_thres), "-k", str(k), "-s", str(s), "-o", msh, skim],
             stderr=open(os.devnull, 'w'))
    return


def temp_func(cov, p, k, l, copy_thres):
    lam = cov * (l - k) / l
    if copy_thres == 1 or p == 1:
        return 1 - np.exp(-lam * p), lam * (1 - p)
    else:
        s = [(lam * p) ** i / np.math.factorial(i) for i in range(copy_thres)]
        return 1 - np.exp(-lam * p) * sum(s), 0


def poisson_survival(jlam, m):
    s = 0
    for t in range(m):
        s += poisson(jlam, t)
    return 1-s


def error_coef(lam, p, m):
    if m == 1:
        return lam * (1 - p)
    else:
        return 0


def binomial_pmf(p, i, j):
    coeff = 1.0 * np.math.factorial(j) / (np.math.factorial(i) * np.math.factorial(j-i))
    return coeff * p**i * (1-p)**(j-i)


def rep_func_2(x, A, B, C, m, lam):
    s = B + C
    for j in range(1, len(A)):
        ss = 0
        for t in range(m):
            for i in range(0, j+1):
              ss += binomial_pmf(x, i, j) * poisson(i*lam, t)
        s += A[j] * ss
    return s


def estimate_dist(sample_1, sample_2, lib_1, lib_2, ce, le, ee, rl, re, k, tran):
    if sample_1 == sample_2 and lib_1 == lib_2:
        return sample_1, sample_2, 0.0
    sample_dir_1 = os.path.join(lib_1, sample_1)
    sample_dir_2 = os.path.join(lib_2, sample_2)
    msh_1 = os.path.join(sample_dir_1, sample_1 + ".msh")
    msh_2 = os.path.join(sample_dir_2, sample_2 + ".msh")
    n_terms = 5
    q_1 = [0]
    q_2 = [0]
    q_1 += re[sample_1]
    q_2 += re[sample_2]
    dist_stderr = check_output(["mash", "dist", msh_1, msh_2], stderr=STDOUT, universal_newlines=True)
    jac = float(dist_stderr.split()[4].split("/")[0]) / float(dist_stderr.split()[4].split("/")[1])
    gl_1 = le[sample_1]
    gl_2 = le[sample_2]
    cov_1 = ce[sample_1]
    cov_2 = ce[sample_2]
    eps_1 = ee[sample_1]
    eps_2 = ee[sample_2]
    p_1 = np.exp(-k * eps_1)
    p_2 = np.exp(-k * eps_2)
    l_1 = rl[sample_1]
    l_2 = rl[sample_2]
    lam_1 = cov_1 * (l_1 - k) / l_1
    lam_2 = cov_2 * (l_2 - k) / l_2
    m_1 = poisson_percent(lam_1 * p_1, 0.95)
    m_2 = poisson_percent(lam_2 * p_2, 0.95)
    r_1 = temp_func(cov_1, p_1, k, l_1, m_1)
    r_2 = temp_func(cov_2, p_2, k, l_2, m_2)
    wp = r_1[0] * r_2[0] * (gl_1 + gl_2) / 2
    zp = sum(r_1) * gl_1 + sum(r_2) * gl_2
    d_0 = max(0, 1 - (zp * jac / (wp * (1 + jac))) ** (1.0 / k))

    A = [0]
    B = 0
    R_1 = 0
    R_2 = 0
    D = 0
    for j in range(1, n_terms + 1):
        A.append(q_1[j] * poisson_survival(j*lam_1*p_1, m_1))
        B -= A[j]
        R_1 += q_1[j] * j
        R_2 += q_2[j] * j
        D += q_1[j] * poisson_survival(j*lam_1*p_1, m_1) + q_2[j] * poisson_survival(j*lam_2*p_2, m_2)

    C = (jac / (1 + jac)) * (D + R_1 * error_coef(lam_1, p_1, m_1) + R_2 * error_coef(lam_2, p_2, m_2))

    x_0 = (1 - d_0) ** k
    x = newton(rep_func_2, x_0, args=(A, B, C, m_2, lam_2*p_2), maxiter=5000)
    if 0 < x < 1:
        d = 1 - x ** (1.0/k)
    elif x > 1:
        d = 0
    else:
        raise ValueError('Distance larger than 1!')

    if tran:
        if d < 0.75:
            d = -0.75 * np.log(1 - 4.0 * d / 3.0)
        else:
            raise ValueError('Distance between {0} and {1} is not in range [0-0.75]; Unable to apply Jukes-Cantor ' +
                     'transformation'.format(sample_1, sample_2))
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
    samples_names = [f.split('.f')[0] for f in files_names]

    # Making a list of genome-skim files
    genome_skims = [os.path.join(args.input_dir, f) for f in files_names]

    # Initializing distance dataframe
    index = pd.MultiIndex.from_product([samples_names, samples_names], names=['sample', 'sample_2'])
    result_df = pd.DataFrame(columns=index)

    # Initializing coverage, genome length, error rate, and read length dictionaries
    cov_est = dict()
    len_est = dict()
    err_est = dict()
    read_len = dict()
    rep_est = dict()

    # Hard-coded param
    coverage_threshold = 5

    # Number of pools and threads for multi-processing
    n_pool = min(args.p, len(genome_skims))
    n_thread_cov = int(args.p / n_pool)
    n_proc_cov = n_pool * n_thread_cov
    n_pool_dist = min(args.p, len(genome_skims) ** 2)

    # Computing coverage, genome length, error rate, and read length
    sys.stderr.write('[skmer] Estimating coverages using {0} processors...\n'.format(n_proc_cov))
    pool_cov = mp.Pool(n_pool)
    results_cov = [pool_cov.apply_async(estimate_cov, args=(gm, args.l, args.k, args.e, n_thread_cov))
                   for gm in genome_skims]
    for result in results_cov:
        (name, coverage, genome_length, error_rate, read_length, repeat_profile) = result.get()
        cov_est[name] = coverage
        len_est[name] = genome_length
        err_est[name] = error_rate
        read_len[name] = read_length
        rep_est[name] = repeat_profile
    pool_cov.close()
    pool_cov.join()

    # Sketching genome-skims
    sys.stderr.write('[skmer] Sketching genome-skims using {0} processors...\n'.format(n_pool))
    pool_sketch = mp.Pool(n_pool)
    results_sketch = [pool_sketch.apply_async(sketch, args=(gm, args.l, cov_est, err_est, args.k, args.s, read_len))
                      for gm in genome_skims]
    for result in results_sketch:
        result.get()
    pool_sketch.close()
    pool_sketch.join()

    # Estimating pair-wise distances
    sys.stderr.write('[skmer] Estimating distances using {0} processors...\n'.format(n_pool_dist))
    pool_dist = mp.Pool(n_pool_dist)
    results_dist = [pool_dist.apply_async(estimate_dist, args=(s1, s2, args.l, args.l, cov_est, len_est, err_est,
                                                               read_len, rep_est, args.k, args.t))
                    for s1 in samples_names for s2 in samples_names]

    for result in results_dist:
        dist_output = result.get()
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
    sample = os.path.basename(args.input).split('.f')[0]
    sample_dir = os.path.join(os.getcwd(), sample)
    try:
        os.makedirs(sample_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    # Making a list of references samples
    refs = [item for item in os.listdir(args.library) if os.path.isdir(os.path.join(args.library, item))]

    # Initializing distances series
    result_s = pd.Series(index=refs, name=sample)

    # Loading coverage, genome length, error rate, and read length information
    cov_est = dict()
    len_est = dict()
    err_est = dict()
    read_len = dict()
    rep_est = dict()
    for ref in refs:
        ref_dir = os.path.join(args.library, ref)
        info_file = os.path.join(ref_dir, ref + '.dat')
        with open(info_file) as f:
            info = f.read()
        cov_est[ref] = float(info.split('\n')[0].split('\t')[1])
        len_est[ref] = float(info.split('\n')[1].split('\t')[1])
        err_est[ref] = float(info.split('\n')[2].split('\t')[1])
        read_len[ref] = float(info.split('\n')[3].split('\t')[1])
        rep_est[ref] = [float(x) for x in info.split('\n')[4].split('\t')[1:]]

    # Hard-coded param
    coverage_threshold = 5

    # Number of pools for multi-processing
    n_pool_dist = min(args.p, len(refs))

    # Computing the coverage, genome length, error rate, and read length of query sample
    sys.stderr.write('[skmer] Estimating the coverage using {0} processors...\n'.format(args.p))
    (dummy, coverage, genome_length, error_rate, read_length, repeat_profile) = estimate_cov(args.input, os.getcwd(), kl,
                                                                                          args.e, args.p)
    cov_est[sample] = coverage
    len_est[sample] = genome_length
    err_est[sample] = error_rate
    read_len[sample] = read_length
    rep_est[sample] = repeat_profile

    # Sketching the query genome-skim
    sys.stderr.write('[skmer] Sketching the genome-skim...\n')
    sketch(args.input, os.getcwd(), cov_est, err_est, kl, ss, read_len)

    # Estimating pair-wise distances
    sys.stderr.write('[skmer] Estimating distances using {0} processors...\n'.format(n_pool_dist))
    pool_dist = mp.Pool(n_pool_dist)
    results_dist = [pool_dist.apply_async(estimate_dist, args=(sample, ref, os.getcwd(), args.library, cov_est, len_est,
                                                               err_est, read_len, rep_est, kl, args.t))
                    for ref in refs]
    for result in results_dist:
        dist_output = result.get()
        result_s[dist_output[1]] = dist_output[2]

    # Writing distances to file
    sys.stderr.write('[skmer] Writing to file...\n')
    result_s.sort_values(inplace=True)
    result_sr = result_s.apply(repr)
    result_sr.to_csv('{0}-{1}.txt'.format(args.o, sample.lower()), sep='\t', mode='w')

    # Adding query to the reference library
    if args.a:
        os.rename(sample_dir, os.path.join(args.library, sample))
    else:
        shutil.rmtree(sample_dir)


def main():
    # Input arguments parser
    parser = argparse.ArgumentParser(description='Estimate gonomic distances between genome-skims',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument('-v', '--version', action='store_true', help='print the current version')
    subparsers = parser.add_subparsers(title='commands',
                                       description='reference   Process a library of reference genome-skims\n'
                                                   'query       Compare an input genome-skim against a reference ' +
                                                   'library',
                                       help='Run skmer {reference,query} [-h] for additional help')

    # Reference command subparser
    parser_ref = subparsers.add_parser('reference', description='Process a library of reference genome-skims')
    parser_ref.add_argument('input_dir', help='Directory of input genome-skims')
    parser_ref.add_argument('-l', default=os.path.join(os.getcwd(), 'library'),
                            help='Directory of output (reference) library. Default: working_directory/library')
    parser_ref.add_argument('-o', default='ref-dist-mat',
                            help='Output (distances) prefix. Default: ref-dist-mat')
    parser_ref.add_argument('-k', type=int, choices=list(range(1, 32)), default=31, help='K-mer length [1-31]. ' +
                                                                                         'Default: 31', metavar='K')
    parser_ref.add_argument('-s', type=int, default=10**7, help='Sketch size. Default: 10000000')
    parser_ref.add_argument('-e', type=float, help='Base error rate. By default, the error rate is automatically '
                                                   'estimated.')
    parser_ref.add_argument('-t', action='store_true',
                            help='Apply Jukes-Cantor transformation to distances')
    parser_ref.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_ref.set_defaults(func=reference)

    # query command subparser
    parser_qry = subparsers.add_parser('query', description='Compare an input genome-skim against a reference library')
    parser_qry.add_argument('input', help='Input (query) genome-skim')
    parser_qry.add_argument('library', help='Directory of (reference) library')
    parser_qry.add_argument('-a', action='store_true',
                            help='Add the processed input (query) to the (reference) library')
    parser_qry.add_argument('-o', default='dist',
                            help='Output (distances) prefix. Default: dist')
    parser_qry.add_argument('-e', type=float, help='Base error rate. By default, the error rate is automatically '
                                                   'estimated.')
    parser_qry.add_argument('-t', action='store_true',
                            help='Apply Jukes-Cantor transformation to distances')
    parser_qry.add_argument('-p', type=int, choices=list(range(1, mp.cpu_count() + 1)), default=mp.cpu_count(),
                            help='Max number of processors to use [1-{0}]. '.format(mp.cpu_count()) +
                                 'Default for this machine: {0}'.format(mp.cpu_count()), metavar='P')
    parser_qry.set_defaults(func=query)

    args = parser.parse_args()

    # if args.version:
    #     print(__version__)

    args.func(args)


if __name__ == "__main__":
    main()
