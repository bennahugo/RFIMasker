#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 SKA South Africa
#
# This file is part of RFIMasker.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from pyrap.tables import table
import os
import argparse
import numpy as np
import logging
import sys

def main():
    logging.basicConfig(format='%(asctime)s:\t%(message)s')
    logging.StreamHandler(sys.stdout)
    log = logging.getLogger("mask_ms")
    log.setLevel(logging.DEBUG)

    log.info("RFI Masker")
    log.info("----------------------------------------")
    parser = argparse.ArgumentParser(description="Masks a measurement set for known RFI contaminated channels")
    parser.add_argument("ms", metavar="ms", type=str, nargs="+",
                        help="specify one or more measurement sets to flag")
    parser.add_argument("-m", "--mask", type=str, required=True,
                        help="a numpy array of shape [channels] containing a boolean "
                             "per channel")
    parser.add_argument("--accumulation_mode", type=str,
                        choices=["or", "override"], default="or",
                        help="specifies whether mask should override current flags or be "
                             "added (or) to the current")

    parser.add_argument("-s", "--statistics", action='store_true',
                        help="Computes and reports some statistics about the flagged RFI in the MS")
    parser.add_argument("--memory", type=int, default=5,
                        help="Maximum memory to consume in MB for the flag buffer")
    args = parser.parse_args()

    # Load mask
    try:
        mask = np.load(args.mask)
    except:
        raise RuntimeError("Mask %s is not a numpy array" % args.mask)
    if mask.dtype != bool:
        raise RuntimeError("Mask is not a boolean array")
    log.info("Mask %s (%d chan) loaded successfully" % (args.mask, len(mask)))

    # Check each measurement set for compatibility before applying
    for ms in args.ms:
        if not os.path.isdir(ms):
            raise RuntimeError("Measurement set %s does not exist. Check input" % ms)
        t = table(ms + "::SPECTRAL_WINDOW", readonly=True, ack=False)
        spw_name = t.getcol("NAME")
        spw_chans = t.getcol("NUM_CHAN")
        t.close()
        for spw, nch in zip(spw_name, spw_chans):
            if nch != len(mask):
                raise RuntimeError("SPW '%s' of %s contains %d channels but the mask has"
                                   " %d channels. Aborting." % (spw, ms, nch, len(mask)))
        t = table(ms, readonly=True, ack=False)
        flag_shape = t.getcell("FLAG", 0).shape

        nrows = t.nrows()
        t.close()
        if len(flag_shape) != 2: #spectral flags are optional in CASA memo 229
            raise RuntimeError("%s does not support storing spectral flags. "
                               "Maybe run pyxis ms.prep?" % ms)
        log.info("%s appears to be a valid measurement set with %d rows" % (ms, nrows))

    # Apply flags
    for ms_i, ms in enumerate(args.ms):
        t = table(ms, readonly=False, ack=False)
        flag_shape = t.getcell("FLAG", 0).shape
        nrows_to_read = int(args.memory / float((flag_shape[0] * flag_shape[1]) * 1.0 / 1024**2))
        nchunk = int(np.ceil(t.nrows() / float(nrows_to_read)))
        preflags_sum = 0
        postflags_sum = 0
        for chunk_i in xrange(nchunk):
            flag_buffer = t.getcol("FLAG",
                                   chunk_i * nrows_to_read,
                                   min(t.nrows() - (chunk_i * nrows_to_read),
                                   nrows_to_read))
            if args.statistics:
                preflags_sum += np.sum(flag_buffer)

            # for now all correlations flagged equal
            mask_corrs = np.repeat(mask,
                                   flag_buffer.shape[2]).reshape([flag_buffer.shape[1],
                                                                  flag_buffer.shape[2]])
            if args.accumulation_mode == "or":
                flag_buffer[:, :, :] |= mask_corrs
            elif args.accumulation_mode == "override":
                flag_buffer[:, :, :] = mask_corrs
            else:
                pass

            if args.statistics:
                postflags_sum += np.sum(flag_buffer)

            t.putcol("FLAG",
                     flag_buffer,
                     chunk_i * nrows_to_read,
                     min(t.nrows() - (chunk_i * nrows_to_read),
                         nrows_to_read))

        if args.statistics:

            pre_flagged = preflags_sum / float(t.nrows() * flag_shape[0] * flag_shape[1]) * 100.0
            post_flagged = postflags_sum / float(t.nrows() * flag_shape[0] * flag_shape[1]) * 100.0
            log.info("[%d / %d]: %s had %f %% flagged visibilities " \
                     "before masking. After masking " \
                     "it has %f %% flagged" % (ms_i + 1,
                                               len(args.ms),
                                               ms,
                                               pre_flagged,
                                               post_flagged))
        else:
            log.info("[%d / %d]: %s has been flagged" % (ms_i + 1,
                                                         len(args.ms),
                                                         ms))
        t.close()
        del flag_buffer
    log.info("RFI Masker terminated successfully")


if __name__ == "__main__":
    main()
