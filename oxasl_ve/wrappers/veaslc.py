"""
Wrapper for epi_reg command
"""

import fsl.utils.assertions as asrt
from fsl.wrappers import wrapperutils  as wutils

@wutils.fileOrImage('data', 'roi', outprefix='out')
@wutils.fileOrArray('veslocs', 'encdef', 'modmat')
@wutils.fileOrText(' ')
@wutils.fslwrapper
def veaslc(data, roi, veslocs, encdef, imlist, modmat, out="veaslc", **kwargs):
    """
    Wrapper for the ``veasl`` command.
    
    Required options:
    
    Additional options:
    """

    valmap = {
        'inferv' : wutils.SHOW_IF_TRUE,
        'debug' : wutils.SHOW_IF_TRUE,
        'diff' : wutils.SHOW_IF_TRUE,
    }

    asrt.assertIsNifti(data)

    cmd = ['veasl', '--data=%s' % data, '--mask=%s' % roi, '--enc-setup=%s' % encdef, 
           '--imlist=%s' % imlist, '--vessels=%s' % veslocs, '--modmat=%s' % modmat,
           '--out=%s' % out]
    if kwargs.pop("method", "map") == "map":
        cmd.append('--map')

    cmd += wutils.applyArgStyle('--=', valmap=valmap, singlechar_args=True, **kwargs)
    return cmd
