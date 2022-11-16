"""
Python interface to VEASLC

Library for processing vessel-encoded ASL data

Copyright (c) 2008-2013 Univerisity of Oxford
"""
import math
import traceback

import numpy as np

from fsl.data.image import Image

from oxasl import basil
from oxasl.options import OptionCategory, IgnorableOptionGroup
from oxasl.reporting import Report

from ._version import __version__
from .veaslc_cli_wrapper import veaslc_wrapper
from .modmat_default import modmat_default

def veslocs_to_enc(veslocs, nvols=8):
    """
    Automatically generate an encoding scheme matrix given initial vessel locations
    """      
    if nvols not in (6, 8):
        raise ValueError("Auto-generation of encoding matrix only supported for 6 or 8 encoding cycles")

    num_vessels = veslocs.shape[-1]
    if num_vessels != 4:
        raise ValueError("Auto-generation of encoding matrix only supported with 4 inferred vessels")

    try:
        lr1, lr2 = veslocs[0, 0], veslocs[0, 1]
        ap1, ap2 = np.mean(veslocs[1, :2]), np.mean(veslocs[1, 2:])
        two = [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
            [90, 2, lr1, lr2],
            [90, 3, lr1, lr2],
            [0, 2, ap1, ap2],
            [0, 3, ap1, ap2],
        ]
        if nvols == 8:
            # Vector from RC to LV
            LV_minus_RC = veslocs[:2, 3] - veslocs[:2, 0]

            # Want to tag RC and LV simultaneously - gradient angle required
            # is acos[normalised(LV - RC).x]
            tag_rad = math.acos(LV_minus_RC[0] / np.linalg.norm(LV_minus_RC))
            tag_deg = math.degrees(tag_rad)

            # Unit vector in gradient direction
            G = [math.sin(tag_rad), math.cos(tag_rad)]

            # Calculate distance from isocentre to each vessel
            # Dot product of location with gradient unit vector
            isodist = [sum(veslocs[:2, v] * G) for v in range(num_vessels)]
            vA = (isodist[0] + isodist[3])/2
            vB = vA + (abs(vA -isodist[1]) + abs(vA - isodist[2]))/2
            two += [
                [tag_deg, 2, vA, vB],
                [tag_deg, 3, vA, vB],
            ]

        return np.array(two, dtype=np.float)
    except: 
        traceback.print_exc()
        raise RuntimeError("Error generating encoding matrix from vessel locations. Check your vessel locations file.")

def two_to_mac(two):
    """
    Convert TWO matrix to MAC format
    """
    # Encoding cycles - second column
    enccyc = two[:, 1] > 1

    # angles TWO measures angle from AP anti-clock, MAC from LR clock
    th = -two[:, 0][enccyc] + 90
    # MAC uses 180 rotation to indicate reversal of modulation function, thus it
    # is important which of vA or Vb is lower, as for TWO this would reverse
    # the modulation function
    th[two[enccyc, 2] > two[enccyc, 3]] = th[two[enccyc, 2] > two[enccyc, 3]] + 180

    # scales
    D = np.abs(two[enccyc, 2] - two[enccyc, 3]) / 2

    # centres
    thtsp = two[enccyc, 0] * 3.14159265 / 180
    conline = np.mean(two[enccyc, 2:4], 1)

    cx = conline * np.sin(thtsp)
    cy = conline * np.cos(thtsp)

    # reverse cycles
    # these are cycles where the tagging of vessels is reversed, Tom does this
    # by shifting the phase of the moudlation function. So this is NOT 180
    # rotation in my convention, but a shift in the encoding centre
    revcyc = two[enccyc, 1] > 2
    cx[revcyc] = (conline[revcyc] + 2*D[revcyc]) * np.sin(thtsp[revcyc])
    cy[revcyc] = (conline[revcyc] + 2*D[revcyc]) * np.cos(thtsp[revcyc])
    
    # imlist
    imlist = two[:, 1] - 1
    inc = 1
    for idx, img in enumerate(imlist[:]):
        if img > 0:
            imlist[idx] = inc
            inc += 1

    mac = [cx, cy, th, D]
    return np.array(mac, dtype=np.float), imlist

def mac_to_two(mac):
    """ 
    Convert MAC format encoding into TWO format

    This is an inverse of two_to_mac done with remarkably little understanding.
    It seems to assume that 'reverse cycles' occur in odd numbered images which 
    seems unreasonable, but I can't see an obvious way to detect this otherwise
    """
    # MAC format does not explicitly describe the tag/control images
    nenc = mac.shape[1] + 2

    th = np.zeros(nenc)
    th_mac = np.zeros(nenc)
    va = np.zeros(nenc)
    vb = np.zeros(nenc)
    cx = np.zeros(nenc)
    cy = np.zeros(nenc)
    d = np.zeros(nenc)
    imlist = np.arange(nenc) - 1
    imtype = np.arange(nenc)
    imtype[imtype > 1] = (imtype[imtype > 1] % 2) + 2

    cx[2:] = mac[0, :]
    cy[2:] = mac[1, :]
    d[2:] = mac[3, :]
    th[2:] = mac[2, :]
    th_mac[2:] = mac[2, :]

    # Angles
    #
    # MAC uses 180 rotation to indicate reversal of modulation function, thus it
    # is important which of vA or Vb is lower, as for TWO this would reverse
    # the modulation function
    # TWO measures angle from AP anti-clock, MAC from LR clock
    rev_mod = th > 180
    th[rev_mod] = th[rev_mod] - 180
    th[2:] = -th[2:] + 90

    # Scales and centres
    for idx in range(nenc):
        if idx > 1:
            s = np.sin(th[idx] * 3.14159265 / 180)
            c = np.cos(th[idx] * 3.14159265 / 180)
            if np.abs(s) > np.abs(c):
                vb[idx] = cx[idx] / s + d[idx]
                va[idx] = cx[idx] / s - d[idx]
            else:
                vb[idx] = cy[idx] / c + d[idx]
                va[idx] = cy[idx] / c - d[idx]

        if idx > 1 and idx % 2 == 1:
            vb[idx] = vb[idx] - 2*d[idx]
            va[idx] = va[idx] - 2*d[idx]
    
        if th_mac[idx] > 180:
            va[idx], vb[idx] = vb[idx], va[idx]

    two = np.column_stack((th, imtype, va, vb))
    return two, imlist

def generate_mask(data, imlist, frac=0.5):
    """
    Auto-generate the mask used for inference of vessel positions/proportions by VEASL

    Note that this is *not* the mask used to output the flows etc

    :param data: ASL data as Numpy array
    :param imlist: Image definition list
    :param frac: Percentile fraction to use

    :return: Numpy integer array defining inference mask
    """
    imlist = list(imlist) # Might by np.ndarray
    tag_idx, ctl_idx = imlist.index(-1), imlist.index(0)
    diffdata = np.abs(data[..., tag_idx] - data[..., ctl_idx])
    thresh = np.percentile(diffdata, 99) * (1-frac)
    return (diffdata > thresh).astype(np.int)

def _autogen_mask(wsp):
    maskdata = generate_mask(wsp.asldata.data, wsp.imlist, wsp.ifnone("infer_mask_frac", 0.5))
    wsp.infer_mask = Image(maskdata, header=wsp.asldata.header)

def _decode_infer(wsp):
    """
    Run VEASLC vessel decoding
    """
    # Generate mask for inference
    if wsp.infer_mask is None:
        _autogen_mask(wsp)

    # Run VEASLC
    flow, prob, extras, log = veaslc_wrapper(wsp, wsp.asldata, wsp.infer_mask)
    wsp.set_item("veasl_log", log, save_fn=str)
    wsp.flow = flow
    wsp.prob = prob
    wsp.pis = np.array(extras["pis"])
    wsp.veslocs = np.array([extras["x"], extras["y"]])
    wsp.log.write("   - Vessel locations (inference: %s):\n" % wsp.infer_loc)
    wsp.log.write("     X: %s\n" % wsp.veslocs[0, :])
    wsp.log.write("     Y: %s\n" % wsp.veslocs[1, :])
    if wsp.infer_loc == "rigid":
        tx, ty, rot = tuple(extras["trans"])
        wsp.log.write("     Translation: %.3g, %.3g  Rotation: %.3g (degrees)\n" % (tx, ty, rot * 180 / math.pi))
    wsp.log.write("   - Class proportions:\n")
    wsp.log.write("     %s\n" % wsp.pis)

def _decode(wsp):
    if wsp.veasl is None:
        wsp.sub("veasl")

    if wsp.veslocs.ndim != 2 or wsp.veslocs.shape[0] != 2:
        raise ValueError("Vessel locations should have two rows (XY co-ordinates) and one column per vessel")
    
    wsp.log.write("\nPerforming vessel decoding\n")
    wsp.log.write(" - Initial vessel locations:\n")
    wsp.log.write("   X: %s\n" % wsp.veslocs[0, :])
    wsp.log.write("   Y: %s\n" % wsp.veslocs[1, :])
    num_vessels = wsp.veslocs.shape[1]

    if wsp.encdef is not None:
        # Encoding definition supplied by user
        wsp.log.write("\n - Encoding definition supplied by user in '%s' format\n" % wsp.encdef_format)
        if wsp.encdef_format == "mac":
            wsp.veasl.enc_mac = wsp.encdef
            wsp.veasl.enc_two, imlist = mac_to_two(wsp.veasl.enc_mac)
        elif wsp.encdef_format == "two":
            wsp.veasl.enc_two = wsp.encdef
            wsp.veasl.enc_mac, imlist = two_to_mac(wsp.veasl.enc_two)
        else:
            raise ValueError("Unknown encoding definition format: %s" % wsp.encdef_format)

        # imlist describes the series of encoded volumes - normally the default is fine
        if wsp.imlist is None:
            wsp.veasl.imlist = imlist
    else:
        # Auto-generate encoding definition from the initial vessel locations
        wsp.log.write("\n - Encoding definition will be automatically generated from vessel locations")
        wsp.veasl.enc_two = veslocs_to_enc(wsp.veslocs, wsp.asldata.nenc)
        wsp.veasl.enc_mac, wsp.veasl.imlist = two_to_mac(wsp.veasl.enc_two)

    wsp.log.write("\n - Encoding matrix:\nTWO\n")
    for row in wsp.veasl.enc_two:
        wsp.log.write("   %s\n" % ", ".join([str(v) for v in row]))
    wsp.log.write("MAC:\n")
    for row in wsp.veasl.enc_mac:
        wsp.log.write("   %s\n" % ", ".join([str(v) for v in row]))

    # Modulation matrix
    if wsp.modmat is None:
        wsp.veasl.modmat = modmat_default

    # Set iaf on the workspace so it is always possible to detect when we are using
    # pairwise subtracted data
    wsp.veasl.iaf = wsp.asldata.iaf

    # Make sure encoding cycles are together in the data and 
    # average over repeats if required FIXME is this wise/required?
    wsp.veasl.asldata_mar = wsp.asldata.mean_across_repeats(diff=False).reorder(out_order="lrt")
   
    wsp.veasl.infer_loc = wsp.ifnone("infer_loc", "rigid")
    if wsp.init_loc and wsp.asldata.ntis > 1:
        wsp.log.write("\n - Doing initial fit for vessel locations using mean data\n")
        asldata_mean = np.zeros(list(wsp.veasl.asldata_mar.data.shape[:3]) + [wsp.asldata.ntc], dtype=np.float)
        for idx in range(wsp.asldata.ntis):
            asldata_mean += wsp.veasl.asldata_mar[..., idx*wsp.asldata.ntc:(idx+1)*wsp.asldata.ntc]
        asldata_mean /= wsp.asldata.ntis
        asldata_mean = Image(asldata_mean, header=wsp.asldata.header)

        wsp_init = wsp.veasl.sub("init")
        wsp_init.asldata = asldata_mean

        _decode_infer(wsp_init)
        wsp.veasl.veslocs_orig = wsp.veasl.veslocs
        wsp.veasl.veslocs = wsp_init.veslocs
        wsp.veasl.infer_loc = wsp.ifnone("infer_loc_pld", "none")

    # Do vessel decoding on each TI/PLD with fixed vessel locations
    for idx in range(wsp.asldata.ntis):
        wsp.log.write("\n - Fitting PLD %i\n" % (idx+1))
        wsp_ti = wsp.veasl.sub("pld%i" % (idx+1))
        wsp_ti.asldata = wsp.veasl.asldata_mar.single_ti(idx)
        _decode_infer(wsp_ti)
        if wsp.asldata.ntis == 1:
            wsp.veasl.veslocs_orig = wsp.veasl.veslocs
            wsp.veasl.veslocs = wsp_ti.veslocs

    wsp.log.write("\nDONE vessel decoding\n")
    return num_vessels

def _model_vessels(wsp, num_vessels):
    from oxasl import oxford_asl
    wsp.log.write("\nProcessing per-vessel decoded images\n\n")
    wsp.sub("output")

    # Generate per-vessel subtracted images. Note that we multiply the flow by 2. Why?
    #
    # "for a standard ASL sequence, in the label condition we measure S - B (where S is 
    #  static tissue and B is blood signal) and in the control condition we measure S + B, 
    #  so control - label gives us 2B, which is what is normally passed into oxford_asl. 
    #  However, when we explicitly define the encoding matrix and do the decoding, we are 
    #  directly estimating B, rather than 2B, so the output of the decoding is half what 
    #  you would expect from a conventional control - label subtraction."
    #    - Thomas Okell
    for vessel in range(num_vessels):
        wsp.log.write("\n - Processing vessel %i\n" % (vessel+1))
        wsp_ves = wsp.veasl.sub("vessel%i" % (vessel+1))
        vessel_data = np.zeros(list(wsp.asldata.shape[:3]) + [wsp.asldata.ntis,], dtype=np.float)
        for ti_idx in range(wsp.asldata.ntis):
            flow = getattr(wsp.veasl, "pld%i" % (ti_idx+1)).flow
            vessel_data[..., ti_idx] = flow.data[..., vessel] * 2
        vessel_img = wsp.veasl.asldata_mar.derived(vessel_data, iaf="diff", order="rt")
        wsp_ves.asldata = vessel_img
        basil.basil(wsp_ves, wsp_ves)
        subname = "vessel%i" % (vessel+1)
        report = Report("Output for vessel %i" % (vessel+1))
        oxford_asl.output_native(wsp.output.sub(subname), wsp_ves, report=report)
        wsp.report.add(subname, report)

def _combine_vessels_sum(wsp, num_vessels, output):
    """
    Generate combined vessel output by summing the contributions
    from individual vessels
    
    Used for perfusion data
    """
    have_output = False
    all_vessel_img = None
    for vessel in range(num_vessels):
        vessel_wsp = getattr(wsp.output, "vessel%i" % (vessel+1))
        vessel_img = getattr(vessel_wsp.native, output, None)
        if vessel_img is not None:
            if all_vessel_img is None:
                all_vessel_img = np.zeros(vessel_img.shape, dtype=np.float)
            all_vessel_img += vessel_img.data
            have_output = True
    if have_output:
        setattr(wsp.output.all_vessels.native, output, Image(all_vessel_img, header=wsp.asldata.header))

def _combine_vessels_weighted(wsp, num_vessels, output, method="weightedperf"):
    """
    Generate combined vessel output by summing the contributions
    from individual vessels weighted by voxelwise perfusion
    
    Used for arrival time and variance/std dev
    """
    vessel_perf = np.zeros(list(wsp.asldata.shape[:3]) + [num_vessels,], dtype=np.float)
    vessel_output = np.zeros(list(wsp.asldata.shape[:3]) + [num_vessels,], dtype=np.float)
    have_output = False
    for vessel in range(num_vessels):
        vessel_wsp = getattr(wsp.output, "vessel%i" % (vessel+1))
        vessel_img = getattr(vessel_wsp.native, output, None)
        if vessel_img is not None:
            vessel_perf[..., vessel] = vessel_wsp.native.perfusion.data
            vessel_output[..., vessel] = vessel_img.data
            have_output = True

    if have_output:
        if method == "singleperf":
            # This selects the arrival time from a single vessel with highest perfusion
            best_vessel = np.argmax(vessel_perf, -1)
            shape3d, t = vessel_output.shape[:-1], vessel_output.shape[-1]
            flat = vessel_output.reshape(-1, t)[np.arange(np.prod(shape3d)), best_vessel.ravel()]
            all_vessel_output = flat.reshape(shape3d)
            wsp.output.all_vessels.native.best_vessel = Image(best_vessel+1, header=wsp.asldata.header)
        elif method == "weightedperf":
            # Tests the arrival time from a weighted average of vessels weighted by perfusion. Note that
            # we need to protect against divide by zero
            total_perf = wsp.output.all_vessels.native.perfusion.data
            total_perf[total_perf == 0] = 1
            all_vessel_output = np.nan_to_num(np.sum(vessel_perf * vessel_output, axis=-1) / total_perf)
        #elif combine_method == "weightedprob":
        #    # Sets the arrival time as a weighted average by vessel probability FIXME not working right now
        #    all_vessel_output = np.sum(prob * vessel_arrival, axis=-1)
        else:
            raise ValueError("Unrecognized combination method for weighted combination: %s" % method)
        setattr(wsp.output.all_vessels.native, output, Image(all_vessel_output, header=wsp.asldata.header))

def _combine_vessels(wsp, num_vessels):
    """
    Generate combined output for all vessels
    """
    from oxasl import oxford_asl
    wsp.log.write("\nGenerating combined images for all vessels\n\n")

    # Generate combined perfusion and aCBV maps over all vessels
    wsp.output.sub("all_vessels")
    wsp.output.all_vessels.sub("native")
    for otype in ("", "_calib", "_std", "_var", "_std_calib", "_var_calib"):
        for oname in ("perfusion", "aCBV", "modelfit", "arrival"):
            output = "%s%s" % (oname, otype)
            if oname != "arrival" and otype  in ("", "_calib"):
                _combine_vessels_sum(wsp, num_vessels, output)                
            else:
                _combine_vessels_weighted(wsp, num_vessels, output, method=wsp.ifnone("arrival_combine", "weightedperf"))

    report = Report("Combined output for all vessels")
    oxford_asl.output_report(wsp.output.all_vessels.native, report=report)
    wsp.report.add("all_vessels", report)

def model_ve(wsp):
    """
    Do vessel decoding and modelling on vessel-encoded ASL data

    :param wsp: Workspace object

    Required workspace attributes
    -----------------------------
      
      - ``asldata`` - ASLImage containing vessel-encoded data
      - ``veslocs`` - Initial vessel locations as 2xN matrix

    Optional workspace attributes
    -----------------------------

      - ``encdef`` - Encoding matrix in MAC format - autogenerated from vessel
                   locations if not supplied
      - ``modmat`` - Modulation matrix as Numpy array FIXME file at the moment

    See ``VeaslOptionGroup`` for other VEASL options

    Workspace attributes updated
    ----------------------------

      - ``veasl``             - Sub-workspace containing VEASL decoding output
      - ``veasl.enc_mac``     - Encoding matrix in MAC format as Numpy array
      - ``veasl.asldata_mar`` - AslImage containing input ASL data averaged across repeats
      - ``veasl.init``        - Sub workspace containing initialization inference output
                                on data averaged over PLDs. Only included when number of PLDs
                                is > 1 and wsp.init_loc is True
      - ``veasl.init.veslocs``  Inferred vessel locations from initialization run
      - ``veasl.pld<n>``      - Sub workspace containing output for PLD number <n>
      - ``veasl.pld<n>.flow`` - Image containing perfusion weighted signal for each vessel
      - ``veasl.pld<n>.prob`` - Image containing probability for each vessel
      - ``veasl.pld<n>.veslocs` - Inferred vessel locations at this PLD
      - ``veasl.vessel<n>``    - Sub workspace containing model fitting output for vessel 
                                 <n> at all PLDs
      - ``output.vessel<n>``   - Sub workspace containing native/structural/standard space 
                                 parameter maps for vessel <n>
      - ``output.all_vessels`` - Sub workspace containing native/structural/standard space 
                                 parameter maps for all vessels combined
    """
    from oxasl import oxford_asl

    num_vessels = _decode(wsp)
    _model_vessels(wsp, num_vessels)
    _combine_vessels(wsp, num_vessels)

    # Re-do registration using all-vessel PWI map.
    oxford_asl.redo_reg(wsp, wsp.output.all_vessels.native.perfusion)

    oxford_asl.output_trans(wsp.output.all_vessels)

    wsp.log.write("\nDONE processing\n")

class VeaslOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing ASL data
    """

    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "oxasl_ve", **kwargs)

    def groups(self, parser):
        ret = []
        g = IgnorableOptionGroup(parser, "VEASL Options")
        g.add_option("--veslocs", help="Vessel locations file", type="matrix")
        g.add_option("--nfpc", help="Number of flows per class", type="int", default=2)
        g.add_option("--veasl-method", help="VEASL inference method", choices=["map", "mcmc"], default="map")
        g.add_option("--infer-v", help="Infer flow velocity", action="store_true", default=False)
        g.add_option("--infer-loc", help="Vessel location inference method to use. For multi-pld data can also infer vessel locations individually using --infer-loc-pld", choices=["none", "xy", "rigid", "affine"], default="rigid")
        g.add_option("--init-loc", help="For multi-PLD data, initialise vessel locations by fitting with mean data", action="store_true", default=False)
        g.add_option("--infer-loc-pld", help="For multi-PLD data when --init-loc is specified, this specifies an optional additional vessel location inference method during individual PLD fitting.", choices=["none", "xy", "rigid", "affine"], default="none")
        g.add_option("--xy-std", help="Prior standard deviation for vessel positions", type="float", default=1.0)
        g.add_option("--rot-std", help="Prior standard deviation for rotation angle (degrees) if using rigid body inference", type="float", default=1.2)
        g.add_option("--v-mean", help="Prior mean flow velocity if using --infer-v", type="float", default=0.3)
        g.add_option("--v-std", help="Prior standard deviation for flow velocity if using --infer-v", type="float", default=0.01)
        g.add_option("--infer-mask-frac", help="Fraction of 99th percentile to use when generating inference mask", type="float", default=0.5)
        g.add_option("--modmat", help="Modulation matrix file")
        g.add_option("--encdef", help="Encoding definition file", type="matrix")
        g.add_option("--encdef-format", help="Encoding definition file format", default="mac")
        g.add_option("--arrival-combine", help="Method for combining arrival time maps", choices=["weightedperf", "singleperf"], default="weightedperf")
        ret.append(g)
        g = IgnorableOptionGroup(parser, "VEASL options for --method=mcmc")
        g.add_option("--njumps", help="Number of parameter jumps", type="int", default=500)
        g.add_option("--burnin", help="Number of jumps before sampling begins", type="int", default=10)
        g.add_option("--sample-every", help="Sampling frequency in jumps", type="int", default=1)
        ret.append(g)
        return ret
