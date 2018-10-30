from oxasl_ve.wrappers import veaslc

def veaslc_wrapper(wsp, data, roi):
    """
    """
    # Run the C code
    log = veaslc(data, roi, 
                 method=wsp.ifnone("veasl_method", "map"),
                 veslocs=wsp.veslocs, 
                 imlist=wsp.imlist, 
                 encdef=wsp.enc_mac,
                 modmat=wsp.modmat, 
                 nfpc=wsp.nfpc, 
                 #infer_loc=wsp.infer_loc, 
                 inferv=wsp.ifnone("infer_v", False), 
                 #xy_std=wsp.ifnone("xy_std", 1), 
                 #v_mean=wsp.ifnone("v_mean", 0.3), 
                 #v_std=wsp.ifnone("v_std", 0.01), 
                 #rot_std=wsp.ifnone("rot_std", 1.2),
                 njumps=wsp.ifnone("num_jumps", 500), 
                 burnin=wsp.ifnone("burnin", 10), 
                 sampleevery=wsp.ifnone("sample_every", 1), 
                 debug=wsp.ifnone("debug", False))

    return ret_flow, ret_prob, extras, log.decode("UTF-8")
