from fsl.wrappers import LOAD
from oxasl_ve.wrappers import veaslc

def imgnums_to_imlist(imgnums):
    """
    Convert image number to imlist string: -1=T, 0-9 as is,
    10+ use letters A, B...

    :param imgnums: Sequence of numbers (int or float type but must be int values)
    :return: String imlist to pass to veaslc
    """
    imlist = ""
    for imgnum in imgnums:
        if imgnum < 0:
            imlist += "T"
        elif imgnum < 10:
            imlist += str(int(imgnum))
        else:
            imlist += chr(ord('A') + int(imgnum) - 10)
    return imlist

def veaslc_wrapper(wsp, data, roi):
    """
    """
    imlist = imgnums_to_imlist(wsp.imlist)
    # Run the C code
    ret = veaslc(data, roi, out=LOAD,
                 diff=wsp.iaf == "vediff",
                 method=wsp.ifnone("veasl_method", "map"),
                 veslocs=wsp.veslocs, 
                 imlist=imlist,
                 encdef=wsp.enc_mac,
                 modmat=wsp.modmat, 
                 nfpc=wsp.nfpc, 
                 inferloc=wsp.infer_loc, 
                 inferv=wsp.ifnone("infer_v", False), 
                 xystd=wsp.ifnone("xy_std", 1), 
                 rotstd=wsp.ifnone("rot_std", 1.2),
                 vmean=wsp.ifnone("v_mean", 0.3), 
                 vstd=wsp.ifnone("v_std", 0.01), 
                 njumps=wsp.ifnone("num_jumps", 500), 
                 burnin=wsp.ifnone("burnin", 10), 
                 sampleevery=wsp.ifnone("sample_every", 1), 
                 debug=wsp.ifnone("debug", False),
                 log=wsp.fsllog)

    flow = ret["out/flow"]
    prob = ret["out/vessel_prob"]
    extras = {
        "pis" : ret["out/pis"], 
        "x" : ret["out/x"], 
        "y" : ret["out/y"],
        "trans" : ret.get("out/trans", None),
    }
    log = ret["out/logfile"]
    return flow, prob, extras, log
