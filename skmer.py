import numpy as np
from scipy.optimize import newton
import argparse
import os
import sys
import errno
import pandas as pd
from subprocess import call, check_output, STDOUT
import multiprocessing as mp
from Bio import SeqIO
import random


def find_len(sps, fasta_dir, fasta_suffix):
    fasta_file = os.path.join(fasta_dir, sps + fasta_suffix)
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq) - record.seq.count("N")
    return sps, total_length


def gen_read(se, sps, gl, reads_dir, l, e, fasta_dir, fasta_suffix, local_dir):
    fasta_file = os.path.join(fasta_dir, sps + fasta_suffix)
    fastq_dir = os.path.join(local_dir, sps)
    try:
        os.makedirs(fastq_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise
    fastq_file = os.path.join(fastq_dir, sps)
    c = 1.0 * se / gl
    q = int(-10 * np.log10(e))
    call(["art_ill umina", "-i", fasta_file, "-l", str(l), "-f", str(c), "-o", fastq_file, "-na",
          "-qL", str(q), "-qU", str(q), "-q"])
    call(["mv", fastq_dir, reads_dir])
    # fastq_file = os.path.join(fastq_dir, sps + "_p_{0}_cov_{1}.fq".format(group, c))
    # fastq_file_dummy = os.path.join(fastq_dir, sps + "_p_{0}_cov_{1}_dummy.fq".format(group, c))
    # n_reads = int(gl * c / l)
    # call(["wgsim", "-1", str(l), "-2", str(l), "-d", "50", "-N", str(n_reads), "-e", str(e),
    #       "-r", "0.0", "-S", str(seed), fasta_file, fastq_file, fastq_file_dummy])
    # call(["rm", fastq_file_dummy])
    

def cov_func(x, r, p, k, l):
    lam = x * (1.0 * (l-k)) / l
    return lam * (p ** 2) * np.exp(-lam * p) - 2 * r * (p * np.exp(-lam * p) + 1 - p)


def estimate_cov(fastq, k, e, temp_dir, ms, local_dir):
    name = fastq.split('/')[-1].split('.f')[0]
    mercnt = os.path.join(local_dir, name + '.jf')
    if ms == 1:
        call(["jellyfish", "count", "-m", str(k), "-s", "100M", "-C", "-o", mercnt, fastq])
    histo_stderr = check_output(["jellyfish", "histo", mercnt], stderr=STDOUT)
    count = [0]
    for item in histo_stderr.split('\n')[:-2]:
        count.append(int(item.split()[1]))
    ind = count.index(max(count[2:]))
    p0 = np.exp(-k * e)
    rl = int(check_output(["awk", 'NR==2{print length($1); exit}', fastq], stderr=STDOUT).split()[0])
    cov_0 = (1.0 / p0) * (1.0 * rl / (rl - k)) * (ind + 1) * count[ind + 1] / count[ind]
    cov = cov_0
    if cov_0 < 2:
        r21 = 1.0 * count[2] / count[1]
        cov = newton(cov_func, cov_0, args=(r21, p0, k, rl))
    elif 4 < cov_0:
        cov = (1.0 / p0) * (1.0 * rl / (rl - k)) * (ind + 2) * count[ind + 2] / count[ind + 1]
    wc_stderr = check_output(["wc", "-l", fastq], stderr=STDOUT)
    tot_seq = rl * int(wc_stderr.split()[0].strip()) / 4
    g_len = tot_seq / cov
    histo_file = os.path.join(temp_dir, name + '.hist')
    with open(histo_file, mode='w') as f:
        f.write(histo_stderr)
    call(["rm", mercnt])
    return name, cov, g_len


def sketch(fastq, ce, k, e, s, temp_dir, fasta_dir, fasta_suffix, cov_thres):
    name = fastq.split('/')[-1].split('.f')[0]
    msh = os.path.join(temp_dir, name)
    cov = ce[name]
    copy_thres = int(cov / cov_thres) + 1
    if cov < cov_thres or e == 0.0:
        call(["mash", "sketch", "-k", str(k), "-s", str(s), "-r", "-o", msh, fastq])
    else:
        call(["mash", "sketch", "-m", str(copy_thres), "-k", str(k), "-s", str(s), "-o", msh, fastq])
    fasta_file = os.path.join(fasta_dir, name + fasta_suffix)
    msh_assembly = os.path.join(temp_dir, name + '_assembly')
    # call(["mash", "sketch", "-k", str(k), "-s", str(s), "-o", msh_assembly, fasta_file])
    # name = fastq.split('/')[-1].split('.f')[0]
    # cov = cov_dict[name]
    # copy_thres = int(cov / cov_thres) + 1
    # name_2 = name + "_2"
    # fastq = fastq.split('.')[0] + "_2.fastq"
    # msh = os.path.join(temp_dir, name_2)
    # if cov < cov_thres or e == 0.0:
    #     call(["mash", "sketch", "-k", str(k), "-s", str(s), "-r", "-o", msh, fastq])
    # else:
    #     call(["mash", "sketch", "-m", str(copy_thres), "-k", str(k), "-s", str(s), "-o", msh, fastq])
    return


def jacc2dist(j, k, gl1, gl2, len_penalty):
    if len_penalty:
        return 1 - ((gl1 + gl2) * j / ((gl1 + gl2) * (1 + j) / 2)) ** (1.0 / k)
    else:
        return 1 - ((gl1 + gl2) * j / (min(gl1, gl2) * (1 + j))) ** (1.0 / k)


def temp_func(cov, p, k, l, cov_thres):
    copy_thres = int(cov / cov_thres) + 1
    lam = cov * (l - k) / l
    if copy_thres == 1 or p == 1:
        return 1 - np.exp(-lam * p), lam * (1 - p)
    else:
        s = [(lam * p) ** i / np.math.factorial(i) for i in range(copy_thres)]
        return 1 - np.exp(-lam * p) * sum(s), 0


def estimate_dist(fastq_1, fastq_2, ce, le, ls, k, e, temp_dir, cov_thres, idx=None, mode=None, assembly_dir=None):
    name_1 = fastq_1.split('/')[-1].split('.f')[0]
    name_2 = fastq_2.split('/')[-1].split('.f')[0]
    if name_1 == name_2:
        return name_1, name_2, idx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    name_assembly_1 = name_1
    name_assembly_2 = name_2
    msh_assembly_1 = os.path.join(temp_dir, name_assembly_1 + "_assembly.msh")
    msh_assembly_2 = os.path.join(temp_dir, name_assembly_2 + "_assembly.msh")
    if mode == 'cov':
        name_assembly_1 = name_1.split('_p_')[0]
        name_assembly_2 = name_2.split('_p_')[0]
        msh_assembly_1 = os.path.join(assembly_dir, name_assembly_1 + "_assembly.msh")
        msh_assembly_2 = os.path.join(assembly_dir, name_assembly_2 + "_assembly.msh")
    dist_stderr = check_output(["mash", "dist", msh_assembly_1, msh_assembly_2], stderr=STDOUT)
    j_assembly = float(dist_stderr.split()[4].split("/")[0]) / float(dist_stderr.split()[4].split("/")[1])
    d_mashp_assembly = float(dist_stderr.split()[2])
    gl_1 = ls[name_assembly_1]
    gl_2 = ls[name_assembly_2]
    d_mashl_assembly = (-1.0 / k) * np.log((gl_1 + gl_2) * j_assembly / (min(gl_1, gl_2) * (1 + j_assembly)))
    d_macp_assembly = jacc2dist(j_assembly, k, gl_1, gl_2, len_penalty=True)
    d_macl_assembly = jacc2dist(j_assembly, k, gl_1, gl_2, len_penalty=False)
    # damjt = damj
    # if 0.0 < damj < 0.75:
    #     damjt = -0.75 * np.log(1 - 4.0 * damj / 3.0)
    #
    # damjtl = damjl
    # if 0.0 < damjl < 0.75:
    #     damjtl = -0.75 * np.log(1 - 4.0 * damjl / 3.0)
    msh_1 = os.path.join(temp_dir, name_1 + ".msh")
    msh_2 = os.path.join(temp_dir, name_2 + ".msh")
    dist_stderr = check_output(["mash", "dist", msh_1, msh_2], stderr=STDOUT)
    j_read = float(dist_stderr.split()[4].split("/")[0]) / float(dist_stderr.split()[4].split("/")[1])
    d_mashp = float(dist_stderr.split()[2])
    egl_1 = le[name_1]
    egl_2 = le[name_2]
    d_mashl = (-1.0 / k) * np.log((egl_1 + egl_2) * j_read / (min(egl_1, egl_2) * (1 + j_read)))
    d_macp = jacc2dist(j_read, k, egl_1, egl_2, len_penalty=True)
    d_macl = jacc2dist(j_read, k, egl_1, egl_2, len_penalty=False)
    cov_1 = ce[name_1]
    cov_2 = ce[name_2]
    p = np.exp(-k * e)
    l_1 = int(check_output(["awk", 'NR==2{print length($1); exit}', fastq_1], stderr=STDOUT).split()[0])
    l_2 = int(check_output(["awk", 'NR==2{print length($1); exit}', fastq_2], stderr=STDOUT).split()[0])
    r_1 = temp_func(cov_1, p, k, l_1, cov_thres)
    r_2 = temp_func(cov_2, p, k, l_2, cov_thres)
    wp = r_1[0] * r_2[0] * (egl_1 + egl_2) / 2
    zp = sum(r_1) * egl_1 + sum(r_2) * egl_2
    d_corrp = 1 - (zp * j_read / (wp * (1 + j_read))) ** (1.0 / k)
    d_corrp_trans = d_corrp
    if 0.0 < d_corrp < 0.75:
        d_corrp_trans = -0.75 * np.log(1 - 4.0 * d_corrp / 3.0)
    wl = r_1[0] * r_2[0] * min(egl_1, egl_2)
    zl = sum(r_1) * egl_1 + sum(r_2) * egl_2
    d_corrl = 1 - (zl * j_read / (wl * (1 + j_read))) ** (1.0 / k)
    d_corrl_trans = d_corrl
    if 0.0 < d_corrl < 0.75:
        d_corrl_trans = -0.75 * np.log(1 - 4.0 * d_corrl / 3.0)
    return name_1, name_2, idx, d_mashp_assembly, d_mashl_assembly, d_macp_assembly, d_macl_assembly, d_mashp, d_mashl,\
        d_macp, d_macl, d_corrp, d_corrp_trans, d_corrl, d_corrl_trans


def estimate_self_dist(fastq_1, ce, le, k, e, temp_dir, cov_thres):
    name_1 = fastq_1.split('/')[-1].split('.f')[0]
    name_2 = name_1 + "_2"
    msh_1 = os.path.join(temp_dir, name_1 + ".msh")
    msh_2 = os.path.join(temp_dir, name_2 + ".msh")
    dist_stderr = check_output(["mash", "dist", msh_1, msh_2], stderr=STDOUT)
    j_read = float(dist_stderr.split()[4].split("/")[0]) / float(dist_stderr.split()[4].split("/")[1])
    dm = float(dist_stderr.split()[2])
    dmj = 1 - (2 * j_read / (1 + j_read)) ** (1.0 / k)
    cov_1 = ce[name_1]
    cov_2 = cov_1
    gl_1 = le[name_1]
    gl_2 = gl_1
    p = np.exp(-k * e)
    l_1 = int(check_output(["awk", 'NR==2{print length($1); exit}', fastq_1], stderr=STDOUT).split()[0])
    l_2 = l_1
    r_1 = temp_func(cov_1, p, l_1, k, cov_thres)
    r_2 = temp_func(cov_2, p, l_2, k, cov_thres)
    w = r_1[0] * r_2[0] * (gl_1 + gl_2) / 2
    z = sum(r_1) * gl_1 + sum(r_2) * gl_2
    dc = 1 - (z * j_read / (w * (1 + j_read))) ** (1.0 / k)
    dt = dc
    if 0.0 < dc < 0.75:
        dt = -0.75 * np.log(1 - 4.0 * dc / 3.0)
    wl = r_1[0] * r_2[0] * min(gl_1, gl_2)
    zl = sum(r_1) * gl_1 + sum(r_2) * gl_2
    dcl = 1 - (zl * j_read / (wl * (1 + j_read))) ** (1.0 / k)
    dtl = dcl
    if 0.0 < dcl < 0.75:
        dtl = -0.75 * np.log(1 - 4.0 * dcl / 3.0)
    return name_1, dm, dmj, dc, dt, dcl, dtl


def save_self_output(prefix, outdir):
    dist_dict = prefix + '_dist_dict'
    s = pd.Series(eval(dist_dict))
    print '{0} distance: max is {1}, min is {2}'.format(prefix, s.max(), s.min())
    mat_path = os.path.join(outdir, prefix + '-dist-dict.csv')
    s.to_csv(mat_path)


def main():
    # Process command line options
    parser = argparse.ArgumentParser(description='Estimate pair-wise distance for low coverage data')
    parser.add_argument('-i', '--input-list', required=True, help='The list of species')
    parser.add_argument('--input-dir', help='The directory of input sequences', required=True)
    parser.add_argument('-f', '--fasta-suffix', help='The extension of FASTA sequences', default='.fa')
    parser.add_argument('-se', '--seq-effort', required=False, type=int, default=10**8, help='The sequencing effort')
    parser.add_argument('-k', '--kmer-size', required=False, type=int, default=31, help='The length of kmers')
    parser.add_argument('-s', '--sketch-size', required=False, type=int, default=10**7, help='The sketch size')
    parser.add_argument('-l', '--read-length', help='Read length', type=int, default=100)
    parser.add_argument('-e', '--base-error', required=False, type=float, default=0.01, help='Base error rate')
    parser.add_argument('--data-dir', help='The directory where the files are stored', required=True)
    parser.add_argument('--local-dir', help='The local directory where the intermediate files are saved', required=True)
    parser.add_argument('-o', '--output-prefix', help='The prefix for output files name', required=True)
    parser.add_argument('-np', '--num-processors', help='Number of processors to use. ' +
                                                        "Default for this machine is %d" % (mp.cpu_count(),),
                        required=False, type=int, default=mp.cpu_count())
    parser.add_argument('-nr', '--num-runs', help='Number of runs', required=False, type=int, default=10)

    args = parser.parse_args()

    # if not os.path.isdir(args.outputDir) or not os.access(args.outputDir, os.W_OK):
    #     sys.exit("Unable to write to output directory %s" % (args.outputDir,))

    if args.num_processors < 1:
        sys.exit('Number of processors to use must be greater than 0')

    # Creating temporary directory for intermediate data
    temp_dir = os.path.join(args.data_dir, "temp")
    try:
        os.makedirs(temp_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    reads_dir = os.path.join(args.data_dir, "reads")
    try:
        os.makedirs(reads_dir)
    except OSError as Error:
        if Error.errno != errno.EEXIST:
            raise

    # Making a list of fastq files
    species_list = []
    with open(args.input_list) as f:
        for line in f:
            species = line.split()[0]
            species_list.append(species)
    species_list = list(set(species_list))
    read_sets = [os.path.join(os.path.join(reads_dir, species), species + '.fq') for species in species_list]
    index = pd.MultiIndex.from_product([species_list, species_list], names=['species1', 'species2'])
    result_df = pd.DataFrame(columns=index)

    # Initializing coverage and genome length dictionaries
    cov_est = dict()
    len_est = dict()
    len_seq = dict()

    # List of distance methods
    dist_methods = ['mashp_assembly', 'mashl_assembly', 'macp_assembly', 'macl_assembly', 'mashp', 'mashl', 'macp',
                    'macl', 'corrp', 'corrp_trans', 'corrl', 'corrl_trans']

    # Hard-coded param
    # FASTAPATH = "/pedigree2/projects/insect_metagenomics/insectref/insectbase"
    # FASTAPATH = "/pedigree2/projects/insect_metagenomics/birdsref/gigadb"
    coverage_threshold = 5

    # Finding the total length of fasta files
    len_switch = 0
    if len_switch == 1:
        sys.stderr.write('[cash] Finding the length of fastas using {0} processors...\n'.format(args.num_processors))
        pool_len = mp.Pool(args.num_processors)
        results_len = [pool_len.apply_async(find_len, args=(species, args.input_dir, args.fasta_suffix)) 
                       for species in species_list]

        for result in results_len:
            (name, genome_length) = result.get()
            len_seq[name] = genome_length

        with open(os.path.join(temp_dir, "fasta_length.txt"), mode='w') as f:
            for name in len_seq.keys():
                f.write('{0}\t{1}\n'.format(name, len_seq[name]))

        pool_len.close()
        pool_len.join()

    # Loading fasta length data
    if len_switch == 0:
        with open(os.path.join(temp_dir, "fasta_length.txt")) as f:
            for line in f:
                len_seq[line.split('\t')[0]] = int(line.split('\t')[1])

    # Generating reads from fasta files
    read_switch = 0
    if read_switch == 1:
        sys.stderr.write('[cash] Simulating reads with ART using {0} processors...\n'.format(args.num_processors))
        pool_read = mp.Pool(args.num_processors)
        results_read = [pool_read.apply_async(gen_read, args=(args.seq_effort, species, len_seq[species], reads_dir,
                                                              args.read_length, args.base_error, args.input_dir,
                                                              args.fasta_suffix, args.local_dir))
                        for species in species_list]
        for result in results_read:
            result.get()
        pool_read.close()
        pool_read.join()

    # Mixing reads
    mix_switch = 0
    if mix_switch == 1:
        sys.stderr.write('[cash] Randomizing coverage...\n')
        dataset = args.data_dir.split('/')[-1].split('-')[0]
        data_parent_dir = '/'.join(args.data_dir.split('/')[:-1])
        data_dir_list = [os.path.join(data_parent_dir, dataset + seq, "reads") for seq in ["-0.1G", "-0.5G", "-1G"]]
        for species in species_list:
            rnd_dir = random.choice(data_dir_list)
            fastq_dir = os.path.join(rnd_dir, species)
            call(["cp", "-r", fastq_dir, reads_dir])

    # Estimating coverage and genome length
    cov_switch = 1
    mer_switch = 1
    if cov_switch == 1:
        sys.stderr.write('[cash] Estimating coverage using {0} processors...\n'.format(args.num_processors))
        pool_cov = mp.Pool(args.num_processors)
        results_cov = [pool_cov.apply_async(estimate_cov, args=(t, args.kmer_size, args.base_error, temp_dir,
                                                                mer_switch, args.local_dir)) for t in read_sets]

        for result in results_cov:
            (name, coverage, genome_length) = result.get()
            cov_est[name] = coverage
            len_est[name] = genome_length

        with open(os.path.join(temp_dir, "cov.txt"), mode='w') as f:
            for name in cov_est.keys():
                f.write('{0}\t{1}\n'.format(name, cov_est[name]))
        with open(os.path.join(temp_dir, "len.txt"), mode='w') as f:
            for name in len_est.keys():
                f.write('{0}\t{1}\n'.format(name, len_est[name]))

        pool_cov.close()
        pool_cov.join()

    # Loading coverage and genome length data
    if cov_switch == 0:
        with open(os.path.join(temp_dir, "cov.txt")) as f:
            for line in f:
                cov_est[line.split('\t')[0]] = float(line.split('\t')[1])
        with open(os.path.join(temp_dir, "len.txt")) as f:
            for line in f:
                len_est[line.split('\t')[0]] = float(line.split('\t')[1])

    # Sketching read sets
    sketch_switch = 1
    if sketch_switch == 1:
        sys.stderr.write('[cash] Sketching read sets using {0} processors...\n'.format(args.num_processors))
        pool_sketch = mp.Pool(args.num_processors)

        results_sketch = [pool_sketch.apply_async(sketch, args=(t, cov_est, args.kmer_size, args.base_error,
                                                                args.sketch_size, temp_dir, args.input_dir,
                                                                args.fasta_suffix, coverage_threshold))
                          for t in read_sets]
        for result in results_sketch:
            result.get()

        pool_sketch.close()
        pool_sketch.join()

    # Estimating pair-wise distances
    dist_switch = 1
    if dist_switch == 1:
        sys.stderr.write('[cash] Estimating distances using {0} processors...\n'.format(args.num_processors))
        pool_dist = mp.Pool(args.num_processors)
        results_dist = [pool_dist.apply_async(estimate_dist, args=(t, tt, cov_est, len_est, len_seq, args.kmer_size,
                                                                   args.base_error, temp_dir, coverage_threshold))
                        for t in read_sets for tt in read_sets]

        for result in results_dist:
            dist_output = result.get()
            result_df[(dist_output[0], dist_output[1])] = [dist_output[3:]]

    # Writing distances to file
    write_switch = 1
    if write_switch == 1:
        sys.stderr.write('[cash] Writing to file...\n')
        result_dfm = pd.melt(result_df)
        result_dfm[dist_methods] = pd.DataFrame(result_dfm['value'].values.tolist())
        result_dfm.drop('value', axis=1, inplace=True)
        result_prefix = args.output_prefix
        result_suffix = ""
        result_dfm.to_csv(result_prefix + result_suffix + ".txt", sep='\t', mode='w')

    # # Estimating self distances
    # self_dist_switch = 0
    # if self_dist_switch == 1:
    #     sys.stderr.write('[cash] Estimating self distances using {0} processors...\n'.format(args.num_processors))
    #     pool_dist = mp.Pool(args.num_processors)
    #
    #     results_dist = [pool_dist.apply_async(estimate_self_dist, args=(t, cov_est, len_est, args.kmer_size,
    #                                                                     args.base_error, temp_dir, coverage_threshold))
    #                     for t in read_sets]
    #
    #     mash_dist_dict = {}
    #     mac_dist_dict = {}
    #     cor_dist_dict = {}
    #     trans_dist_dict = {}
    #     corl_dist_dict = {}
    #     transl_dist_dict = {}
    #
    #     for result in results_dist:
    #         (name_1, mash_dist, mac_dist, cor_dist, trans_dist, corl_dist, transl_dist) = result.get()
    #         mash_dist_dict[name_1] = mash_dist
    #         mac_dist_dict[name_1] = mac_dist
    #         cor_dist_dict[name_1] = cor_dist
    #         trans_dist_dict[name_1] = trans_dist
    #         corl_dist_dict[name_1] = corl_dist
    #         transl_dist_dict[name_1] = transl_dist
    #
    #     for method in ['mash', 'mac', 'cor', 'trans', 'corl', 'transl']:
    #         save_self_output(method, temp_dir)


if __name__ == "__main__":
    main()
