#!/usr/bin/env python
###############################################################################
#                                                                             #
#    symCalc.py                                                               #
#                                                                             #
#    Toolbox for calculating number of reads to extract from files when       #
#    making simulated datasets                                                #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Release"

###############################################################################

from biom import (parse, table)
from random import random
import json

from parsem import ContigParser

###############################################################################

class BadTotalLengthError(BaseException): pass
class BadreferenceLengthError(BaseException): pass
class BadNumSamplesError(BaseException): pass
class BadNumReferencesError(BaseException): pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Calculator():
    """Everything happens here"""
    def __init__(self,
                 otu,
                 volumes):
        self.bt = self.loadOTUTable(otu)    # biom formatted OTU table
        self.refLengths = []                # lengths of each reference
        self.refUrls = []                   # url (dentifier) for each reference
        self.numRefs = 0                    # len(self.refUrls)
        self.sampleVolumes = volumes        # TOTAL number of reads for each sample
        self.numSamples = len(volumes)

    def calc(self):
        """Parse input files and work out read counts

        returns: A json formatted string containing reference urls and associated counts
        """

        # parse through the OTU table and work out the relative abundances of each ref in each sample
        otu_subset = []     # otu_subset[X] = [ref1 count, ref2 count, ...] for sample X
        done = 0
        for samp_v, samp_id, samp_md in self.bt.iterSamples():
            try:
                sample = samp_v[:self.numRefs]
            except IndexError:
                raise BadNumReferencesError("You want to generate counts for %d references but the OTU table '%s' only contains %d OTUs" % (self.numRefs, otu, len(samp_v)))

            totes = sum(sample)
            otu_subset.append([float(i)/float(totes) for i in sample])

            done += 1
            if done == self.numSamples:
                break

        if done < self.numSamples:
            raise BadNumSamplesError("You want to generate counts for %d samples but the OTU table '%s' only contains %d sites" % (self.numSamples, otu, done))

        # correct for sequence length
        min_bases = float(min(self.refLengths))
        if min_bases == 0:
            raise BadreferenceLengthError("Sample length of 0 is not allowed. Please revise counts in %s" % fileInfo)
        multiplier = [float(i) / min_bases for i in self.refLengths]

        # work out the number to cut
        cuts = []
        for i in range(self.numSamples):
            combined_ratios = [float(multiplier[j] * otu_subset[i][j]) for j in range(self.numRefs)]
            sum_combined = float(sum(combined_ratios))
            final_ratios = [j/sum_combined for j in combined_ratios]
            cuts.append(self.makeReadCounts(final_ratios, self.sampleVolumes[i]))

        ret_dict = []
        for i in range(self.numSamples):
            tmp = {}
            for j in range(self.numRefs):
                tmp[self.refUrls[j]] = cuts[i][j]
            ret_dict.append(tmp)

        return json.dumps(ret_dict)

    def makeReadCounts(self, ratios, volume):
        """split the 'volume' amongst several groups according to the given ratios

        assumes that ratios sums to one
        """
        # for the greater part we can assign reads according to the ratios
        ret_counts = [int(i * volume) for i in ratios]

        # the remainder can get done dynamically
        rem_volume = int(volume - sum(ret_counts))
        # now we assign reads randomly
        if rem_volume > 0:
            # http://stackoverflow.com/questions/3407414/python-sort-without-lambda-expressions
            sorted_indices = [i for (v, i) in sorted((v, i) for (i, v) in enumerate(ratios))]
            sorted_ratios = [ratios[i] for  i in sorted_indices]
            for i in range(1,self.numRefs):
                sorted_ratios[i] += sorted_ratios[i-1]
            # make sure it's one!
            sorted_ratios[-1] = 1.

            for i in range(rem_volume):
                r = random()
                for j in range(self.numRefs):
                    if r <= sorted_ratios[j]:
                        ret_counts[sorted_indices[j]] += 1
                        break

        return ret_counts

    def loadFileInfo(self, fileName):
        """Load information about the reference files

        tab separated

        url1    total_length1
        url2    total_length2
        url3    total_length3
        ...
        """
        lengths_included = False
        with open(fileName, 'r') as fh:
            lc = 0
            for line in fh:
                lc += 1
                if line[0] != '#':
                    fields = line.rstrip().split('\t')
                    self.refUrls.append(fields[0])
                    if len(fields) > 1:
                        try:
                            self.refLengths.append(int(fields[1]))
                        except ValueError:
                            raise BadTotalLengthError("Total lengths must be integers. '%s' on line %d in %s is not valid" % (fields[1], lc, fileName))
                        lengths_included = True

        if not lengths_included:
            # user did not supply contig lengths
            CP = ContigParser.ContigParser()
            for ref in self.refUrls:
                with open(ref, 'r') as cfh:
                    lens = 0
                    for header, seq in CP.readFasta(cfh):
                        lens += len(seq)
                    self.refLengths.append(lens)

        self.numRefs = len(self.refUrls)

    def loadOTUTable(self, fileName):
        """parse the given otu table"""
        raw_file = open(fileName, 'r')
        try:
            # first see if it is in biom format
            return parse.parse_biom_table(raw_file)
        except ValueError:
            # nope, let's hope it's a valid otu table
            raw_file.close()
            raw_file = open(fileName, 'r')
            lines = []
            for line in raw_file:
                lines.append(line)

            sample_ids, obs_ids, data, t_md, t_md_name = parse.parse_classic_table(lines)
            return parse.table_factory(data, sample_ids, obs_ids, None, None)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
