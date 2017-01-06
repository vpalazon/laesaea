# -*- coding: utf-8 -*-

# ===================================================================================

#   File: ba_pivots_BBEd_k.py

# ===================================================================================

#     Copyright (C) 2016 Vicente Palazón-González

#     This program is free software; you can redistribute it and/or modify it under the
#     terms of the GNU General Public License as published by the Free Software
#     Foundation; either version 3 of the License, or (at your option) any later
#     version.

#     This program is distributed in the hope that it will be useful, but WITHOUT ANY
#     WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#     PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License along with this
#     program (file LICENSE.txt); if not, write to the Free Software Foundation, Inc.,
#     675 Mass Ave, Cambridge, MA 02139, USA.

#     You can contact the author at:

#     palazon@uji.es
    
# ===================================================================================

# pivots

import os, sys
from os.path import join

home = os.environ["HOME"]
sys.path[:0] = [join(home, 'eva/scripts')]
from make import system, system_p, make, makeDir, makeF, p

from socket import gethostname
import getopt
hostname = gethostname()

INF = float('infinity')

if __name__ == '__main__':

    def usage():
        sys.stderr.write('use: %s corpus [list of sizes]\n')
        sys.stderr.write('     e.g. corpus = nist8 or random8\n')
        sys.stderr.write('     e.g. [list of sizes] = nist8 nist81000 nist82000...\n')
        sys.exit(1)

    if len(sys.argv) == 1:
        usage()

    pivots_dir = "kt_pivots"

    try:
        os.mkdir(pivots_dir)
    except:
        pass
    
    corpus_name = sys.argv[1]
    lcorpus = []
    for e in sys.argv[2:]:
        lcorpus.append(e)

    if lcorpus == []:
        usage()
        
    lsamples = ['br']

    lpoints = [32]

    lsampling = ['s-sr']

    ## lfunctions = [[function, data_type, vector_size, local_distance]]
    lfunctions = [
        ['f-cc', 'floats', 1, 'levenshtein'],
        ]

    ## lmethods = [[experiment, subexperiment, method, k, step, pivots_start, pivots_end, pivots_step], ...]]
    lexperiment = []
    ## for k in [1, 5, 10, 20, 40]:
    for k in [1]:
        ## lexperiment.append(['klaesa', 'klaesa1', 'BBEd', k, 4, 5, 400, 10])
        lexperiment.append(['klaesa', 'klaesa1', 'BBEd', k, 1, 5, 100, 5])

    for samples in lsamples:
        for corpus in lcorpus:
            for sampling in lsampling:
                for points in lpoints:
                    fn = corpus+samples
                    fn += '_%s%03d' % (sampling, points)
                    for f in lfunctions:
                        function = f[0]
                        data_type = f[1]
                        vector_size = f[2]
                        local_distance = f[3]

                        ext = data_type

                        for experiment in lexperiment:

                            exp = experiment[0]
                            subexp = experiment[1]
                            method = experiment[2]
                            k = experiment[3]
                            step = experiment[4]
                            pivots_start = experiment[5]
                            pivots_end = experiment[6]
                            pivots_step = experiment[7]

                            fn_tt_bin_wext = '%s_%s_%s_%s_%s' % (fn, function, ext, method, local_distance)

                            l_pbc = []

                            for pivots in xrange(pivots_start, pivots_end, pivots_step):
                                fn_kt_wext = fn_tt_bin_wext + '_%s_k%02d_p%02d_i%04d_%s' % (subexp, k, step, pivots, hostname)

                                system_p("./tc -e %s -u %s -m %s -l %s -s %d -d %s -p %d -k %d -f %s_%s.%s -t %s.tt -i %d > kt_pivots/%s.kt" % (exp, subexp, method, local_distance, vector_size, data_type, step, k, fn, function, ext, fn_tt_bin_wext, pivots, fn_kt_wext))
                                for line in file("%s/%s.kt" % (pivots_dir, fn_kt_wext)):
                                    if line.startswith("bound_count:"):
                                        bc = int(line.strip().split()[1].split('/')[0])
                                        break
                                l_pbc.append([pivots, bc])


                            l_pbc.sort()

                            fn_pivots_wext = fn_tt_bin_wext + '_%s_k%02d_p%02d_%s' % (subexp, k, step, hostname)
                            fn_pivots = open(join("./", "%s.pivots" % fn_pivots_wext), "w")
                            for e in l_pbc:
                                print >>fn_pivots, e[0], e[1]
                            fn_pivots.close()

    p('bye.')
