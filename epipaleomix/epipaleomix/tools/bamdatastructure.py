from epipaleomix.tools.commonutils import check_path
import os


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


class MakeCollect(object):
    def __init__(self, make):
        self.makefile = make.pop('Makefile', {})
        self.prefix = self.makefile.pop('Prefixes', {})
        self.beddata = self.makefile.pop('BedFiles', {})
        self.bedfiles, self.bed_plot = self._splitbedopts()
        
    def _splitbedopts(self):
        ''' Split bedfiles options between bedname path and plot boolean '''
        bedf, bedp = {}, {}
        mappacorr = True if self.beddata['EnabledFilter'] else False
        for bedn, bedopts in self.beddata.iteritems():
            if isinstance(bedopts, dict):
                if mappacorr:  ## this is for the outputname
                    bedn += 'MappaOnly'
                bedf[bedn] = bedopts["Path"]
                bedp[bedn] = bedopts["MakeMergePlot"]
            else:  # just options to the specific bedfile
                bedf[bedn] = bedopts
        return bedf, bedp

    
class BamCollect(object):
    def __init__(self, config, bam_name, opts, d_make):
        self.bam_name = bam_name
        self.i_path = os.path.join(config.temp_root, self.bam_name)
        self.o_path = os.path.join(config.destination, self.bam_name)
        check_path(self.i_path)
        check_path(self.o_path)
        self.opts = opts
        self.baminfo = self.opts['BamInfo']
        self.fmt = '{}_{}_{}.txt.gz'
        self.prefix = d_make.prefix  # comes from MakeCollect class

    def retrievedat(self, anal):
        return self._coerce_to_dic(self.opts[anal], self.baminfo,
                                   self.prefix)

    def _coerce_to_dic(self, *args):
        dic = {}
        for arg in args:
            if not isinstance(arg, dict):
                arg = dict([arg])
            dic.update(arg)
        return dic
