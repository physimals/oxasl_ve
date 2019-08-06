from fsl.wrappers import LOAD
from oxasl_ve.wrappers import veaslc

def veaslc_wrapper(wsp, data, roi):
    """
    """
    # Run the C code
    ret = veaslc(data, roi, out=LOAD,
                 diff=data.iaf == "vediff",
                 method=wsp.ifnone("veasl_method", "map"),
                 veslocs=wsp.veslocs, 
                 imlist="T0123456", # FIXME
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
