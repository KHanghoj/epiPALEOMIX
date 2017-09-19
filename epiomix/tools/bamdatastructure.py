from __future__ import print_function
from epiomix.tools.commonutils import check_path
import os


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


class MakeCollect(object):
    def __init__(self, make):
        self.makefile = make.get('Makefile', {})
        self.prefix = self.makefile.get('Prefixes', {})
        self.beddata = self.makefile.get('BedFiles', {})
        self._add_bedname_suffix()

    def _add_bedname_suffix(self):
        ''' Split bedfiles options between bedname path and plot boolean '''
        self.bedfiles = {}
        for bedn, bedopts in self.beddata.items():
            if (isinstance(bedopts, str) and bedopts.endswith(".bed") and
                  self.beddata.get('MappabilityFilter', False)):
                bedn += 'MappaOnly'
            self.bedfiles[bedn] = bedopts

                
    def getfilterinfo(self):
        enabl_filter = self.bedfiles.get('MappabilityFilter', False)
        uniqueness = self.bedfiles.get('MappabilityScore')
        mappapath = self.prefix.get('--MappabilityPath')
        if enabl_filter is True:
            try:
                assert os.path.exists(mappapath), \
                    ("the --MappabilityPath provided '%s' does not exists"
                     % mappapath)
            except TypeError:  # If Mappapath is None 
                raise MakefileError("No '--MappabilityPath' is provided."
                                    " Either add a valid path or disable "
                                    "'MappabilityFilter' and 'GCcorrect'")
        return (enabl_filter, uniqueness, mappapath)


class BamCollect(object):
    def __init__(self, config, bam_name, opts, d_make):
        self.bam_name = bam_name
        
        self.bam_temp_local = os.path.join(config.temp_local, self.bam_name)
        self.bam_output = os.path.join(config.makefiledest, self.bam_name)
        check_path(self.bam_temp_local)
        check_path(self.bam_output)
        self.opts = opts
        self.baminfo = self.opts['BamInfo']
        self.fmt = '{}_{}_{}.txt.gz'
        self.prefix = d_make.prefix  # comes from MakeCollect class
        self._checkindex(self.baminfo['BamPath'], '.bai',
                         'Bam needs to index.\n\t samtools index {}\n')
        self._checkindex(self.prefix['--FastaPath'], '.fai',
                         'Reference needs to index.\n\t samtools faidx {}\n')

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

    def _checkindex(self, path, ext, errormsg):
        test1 = os.path.exists(path+ext)
        bname, _ = os.path.splitext(path)
        test2 = os.path.exists(bname+ext)
        assert (test1+test2) > 0, errormsg.format(path)
        # assert os.path.exists(path+ext), \
        #     errormsg.format(path)
