def clustalo(inf, outdir="default"):
    outfile = inf.replace(".fasta", "_clustalo.aln")
    import os
    if outdir == "default":
        my_clustalo_cmd = "clustalo -i " + inf + " -o " + outfile + " --outfmt clustal"
    else:
        my_clustalo_cmd = "clustalo -i " + inf + " -o " + outdir + "/" + outfile + " --outfmt clustal"
    print(my_clustalo_cmd)
    os.system(my_clustalo_cmd)
    return outfile


def t_coffee(inf, outdir="default"):
    outfile = inf.replace(".fasta", "_tcoffee.aln")
    if outdir == "default":
        from Bio.Align.Applications import TCoffeeCommandline
        tcoffee_cline = TCoffeeCommandline(infile=inf,
                                           output="clustalw",
                                           outfile=outfile)
    else:
        from Bio.Align.Applications import TCoffeeCommandline
        tcoffee_cline = TCoffeeCommandline(infile=inf,
                                           output="clustalw",
                                           outfile=outdir + outfile)
    tcoffee_cline()
    return outfile


def mafft(inf, outdir="default"):
    from Bio.Align.Applications import MafftCommandline
    from io import StringIO
    from Bio import AlignIO
    mafft_cline = MafftCommandline("mafft", input=inf)
    print(mafft_cline)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    outfile = inf.replace('.fasta', '_mafft.aln')
    if outdir == "default":
        AlignIO.write(align, outfile, "clustal")
    else:
        AlignIO.write(align, outdir + '/' + outfile, "clustal")
    return outfile


def clustalw(inf):
    import os
    clustalw_exe = "/Users/tianyilu/Tools/MSA/clustalw-2.1-macosx/clustalw2"
    outfile = inf.replace('.fasta', '_clustalw.aln')
    my_clustal_cmd = clustalw_exe + " -INFILE=" + inf + " -OUTFILE=" + outfile
    os.system(my_clustal_cmd)


def pbcn(inf):
    from Bio.Align.Applications import ProbconsCommandline
    pbcn_exe = "/Users/tianyilu/Tools/MSA/probcons/probcons"
    probcons_cline = ProbconsCommandline(pbcn_exe, input=inf, clustalw=True)
    print(probcons_cline)
    stdout, stderr = probcons_cline()
    with open(inf.replace('.fasta', '_probcons.aln'), "w") as handle:
        handle.write(stdout)
