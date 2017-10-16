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
import re
from RFIMasker import version
from scipy.ndimage.morphology import binary_dilation

pckgdir = os.path.dirname(os.path.abspath(__file__))

def main():
    logging.basicConfig(format='%(asctime)s:\t%(message)s')
    logging.StreamHandler(sys.stdout)
    log = logging.getLogger("mask_ms")
    log.setLevel(logging.DEBUG)

    log.info("RFI Masker")
    log.info("Version %s" % version.__version__)
    log.info("Module installed at %s" % pckgdir)
    log.info("----------------------------------------")
    parser = argparse.ArgumentParser(description="Masks a measurement set for known RFI contaminated channels (version %s)" % (version.__version__))
    parser.add_argument("ms", metavar="ms", type=str, nargs="+",
                        help="specify one or more measurement sets to flag")
    parser.add_argument("-m", "--mask", type=str, required=True,
                        help="a numpy array of shape [channels] containing a boolean "
                             "per channel")
    parser.add_argument("--accumulation_mode", type=str,
                        choices=["or", "override"], default="or",
                        help="specifies whether mask should override current flags or be "
                             "added (or) to the current")
    parser.add_argument("--spwid", nargs="+", type=int, default=[0],
                        help="SPW id (or ids if multiple MSs have been specified")
    parser.add_argument("--dilate", type=str, default=None,
                        help="Dilate mask. This will extend masked regions by this width."
                             "The width can be specified as a frequency band (e.g 0.2MHz),"
                             "or as the the number of channels (20) to dilate by.")

    def _casa_style_range(val):
        """ returns None or tupple with lower and upper bound """
        if not isinstance(val, str):
            raise argparse.ArgumentTypeError("Value must be of type string")
        if val == "":
            return None
        elif re.match(r"^(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?~(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$", val):
            return map(float,val.split("~"))
        else:
            raise argparse.ArgumentError("Value must be range or blank")

    parser.add_argument("--uvrange", default="", type=_casa_style_range,
                        help="UV range to select for flagging in meters (format: lower~upper). Leave blank to select entire array.") 
    parser.add_argument("-s", "--statistics", action='store_true',
                        help="Computes and reports some statistics about the flagged RFI in the MS")
    parser.add_argument("--memory", type=int, default=512,
                        help="Maximum memory to consume in MB for the flag buffer")
    parser.add_argument("--simulate", action='store_true',
                        help="Simulate only, don't store flags")
    args = parser.parse_args()

    # Load mask
    try:
        mask = np.load(args.mask)
        if mask.dtype[0] != np.bool or \
           mask.dtype[1] != np.float64:
           raise RuntimeError("Invalid")
    except:
        raise RuntimeError("Mask %s is not a valid static mask with labelled channel axis [dtype == (bool, float64)]" % args.mask)
    log.info("Mask %s (%d chan) loaded successfully" % (args.mask, len(mask)))
    mask_chans = mask["chans"][1]
    mask_flags = mask["mask"][0]
    # Dilate mask
    if args.dilate:
        try:
            dilate_width = int(args.dilate)
        except ValueError:
            value,units = re.match(r"([\d.]+)([a-zA-Z]+)", args.dilate, re.I).groups()
            if units == 'GHz':
                value = float(value)*1e9
            elif units == 'MHz':
                value = float(value)*1e6
            elif units == 'kHz':
                value = float(value)*1e3
            elif units == 'Hz':
                value = float(value)
            else:
                raise RuntimeError('Unrecognised units for --dilate value::  %s'%units)
        
            chan_width = mask_chans[1] - mask_chans[0]
            dilate_width  = int(value/chan_width) + 1 
        log.info('Dilating mask. Will Flag %d channels on either side of every masked region'%(dilate_width))
        dstruct = np.array([True,True,True])
        mask_flags = binary_dilation(mask_flags, dstruct, iterations=dilate_width)

    masked_channels = mask_chans[np.argwhere(mask_flags)]
    log.info("Mask contains %d masked channels (%.2f %%)" % (len(masked_channels), len(masked_channels)*100.0/float(len(mask_flags))))

    # Check each measurement set for compatibility before applying
    app_spws = args.spwid
    nms = len(args.ms)
    if len(app_spws) == 1:
        app_spws = app_spws * nms
    elif len(app_spws) != nms:
        raise ValueError("You have multiple measurement sets, but specified SPW ids only for some of them. Check input")
    msants = {}
    msddidsel = []
    ms_masks = []
    for msspw, ms in zip(app_spws, args.ms):
        if not os.path.isdir(ms):
            raise RuntimeError("Measurement set %s does not exist. Check input" % ms)

        with table(ms + "::SPECTRAL_WINDOW", readonly=True, ack=False) as t:
            spw_name = t.getcell("NAME", msspw)
            spw_chans = t.getcell("NUM_CHAN", msspw)
            spw_chanlabels = t.getcell("CHAN_FREQ", msspw)
            spw_chanwidths = t.getcell("CHAN_WIDTH", msspw)
            spw_chanlb = spw_chanlabels - spw_chanwidths * 0.5
            spw_chanub = spw_chanlabels + spw_chanwidths * 0.5
            # compute overlap of mask and ms spw channels
            ms_mask = np.sum(np.logical_and(masked_channels > spw_chanlb,
                                            masked_channels < spw_chanub),
                             axis=0) > 0
            ms_masks.append(ms_mask)
            num_mschansmasked = len(np.argwhere(ms_mask))
            log.info("MS '%s' SPW %d contains %d channels falling within the static masked regions (%.2f %%)" %
                     (ms, msspw, num_mschansmasked, num_mschansmasked * 100.0 / spw_chans))

        with table(ms + "::DATA_DESCRIPTION", readonly=True, ack=False) as t:
            msddidsel.append(np.argwhere(t.getcol("SPECTRAL_WINDOW_ID") == msspw)[0])

        with table(ms, readonly=True, ack=False) as t:
            flag_shape = t.getcell("FLAG", 0).shape
            nrows = t.nrows()

        with table(ms + "::ANTENNA", readonly=True, ack=False) as t:
            msants[ms] = t.getcol("POSITION")

        if len(flag_shape) != 2: #spectral flags are optional in CASA memo 229
            raise RuntimeError("%s does not support storing spectral flags. "
                               "Maybe run pyxis ms.prep?" % ms)
        log.info("%s appears to be a valid measurement set with %d rows" % (ms, nrows))

    # Apply flags
    for ms_i, ms in enumerate(args.ms):
        log.info("Processing '%s' (%d / %d)" % (ms, ms_i + 1, len(args.ms)))
        selspw = msddidsel[ms_i]
        t = table(ms, readonly=False, ack=False)
        flag_shape = t.getcell("FLAG", 0).shape
        row_size = flag_shape[0] * flag_shape[1] + 8 * 3 + 8
        nrows_to_read = int(args.memory / float( row_size / 1024.0**2))
        nchunk = int(np.ceil(t.nrows() / float(nrows_to_read)))
        preflags_sum = 0
        postflags_sum = 0
        for chunk_i in xrange(nchunk):
            log.info("\tFlagging chunk %d / %d" % (chunk_i + 1, nchunk))
            nrows_chunk = min(t.nrows() - (chunk_i * nrows_to_read),
                          nrows_to_read)

            flag_buffer = t.getcol("FLAG",
                                   chunk_i * nrows_to_read,
                                   nrows_chunk)
            a1 = t.getcol("ANTENNA1",
                          chunk_i * nrows_to_read,
                          nrows_chunk)
            a2 = t.getcol("ANTENNA2",
                          chunk_i * nrows_to_read,
                          nrows_chunk)
            ddid = t.getcol("DATA_DESC_ID",
                            chunk_i * nrows_to_read,
                            nrows_chunk)
            d2 = np.sum((msants[ms][a1][:] - msants[ms][a2][:])**2, axis=1)

            # ECEF antenna coordinates are in meters. The transforms to get it into UV space are just rotations
            # can just take the euclidian norm here - optimized by not doing sqrt
            luvrange = min(args.uvrange[0], args.uvrange[1]) if args.uvrange is not None else 0.0
            uuvrange = max(args.uvrange[0], args.uvrange[1]) if args.uvrange is not None else np.inf
            sel = np.argwhere(np.logical_and(ddid == selspw,
                                             np.logical_and(d2 >= luvrange**2,
                                                            d2 <= uuvrange**2)))
            log.info("\t\tCriteria resulted in %d out of %d rows being selected for this chunk" % (len(sel),nrows_chunk))
            if args.statistics:
                preflags_sum += np.sum(flag_buffer)

            # for now all correlations flagged equal
            mask_corrs = np.repeat(ms_masks[ms_i],
                                   flag_buffer.shape[2]).reshape([flag_buffer.shape[1],
                                                                  flag_buffer.shape[2]])
            if args.accumulation_mode == "or":
                flag_buffer[sel, :, :] |= mask_corrs
            elif args.accumulation_mode == "override":
                flag_buffer[sel, :, :] = mask_corrs
            else:
                pass

            if args.statistics:
                postflags_sum += np.sum(flag_buffer)
            if not args.simulate:
                log.info("\t\tSaving chunk...")
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
