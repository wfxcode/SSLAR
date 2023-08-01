from STRIDE.Deconvolution import Deconvolve

def R_STRIDE(sc_count_file,sc_anno_file,st_count_file,out_dir,out_prefix):

    normalize = True

    res = Deconvolve(sc_count_file=sc_count_file, sc_anno_file=sc_anno_file, st_count_file=st_count_file,
                     out_dir=out_dir, out_prefix=out_prefix, normalize=normalize,
                     sc_scale_factor=None, st_scale_factor=None)
    return res
