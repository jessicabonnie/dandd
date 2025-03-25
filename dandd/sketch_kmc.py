from __future__ import annotations
import os
import subprocess
import glob
from dandd.sketch_base import SketchObj
from dandd.utils import canon_command
import tempfile



class KMCSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches)
    
    def parse_card(self, proc):
        '''Parse the cardinality streaming from standard out for the kmc cardinality command (kmc_tools info)'''
        for line in proc.stdout.splitlines():
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            
            if ',' in line:  # Handle comma-separated format from card_command
                key, value = line.split(',')
                self.speciesinfo.cardkey[key.strip()] = float(value.strip())
            
            elif 'total k-mers' in line:  # Handle formatted output from kmc_tools info
                try:
                    value = line.split(':')[1].strip()
                    self.speciesinfo.cardkey[self.sfp.full] = float(value)
                except (IndexError, ValueError) as e:
                    print(f"Error parsing line: {line}")
                    raise

    def card_command(self, sketch_paths:list=[]) -> str:
        '''Create the bash command to capture the cardinality of the databases'''
        if len(sketch_paths) == 0:
            if self.kval == 0:
                return
            sketch_paths=[self.sketch]
        cmd= "for db in "+ " ".join(sketch_paths) + "; do value=$(kmc_tools -hp info $db | grep 'total k-mers' | sed 's/ //g' | sed 's/totalk-mers://g'); echo $db,$value; done"
        #kmc_tools info $sketchdir/ngen1/k10/allvar_HG00171_1.fasta.gz_k10 | grep 'total k-mers' | sed 's/ //g' | sed "s/totalk-mers://
        # cmdlist = ["kmc_tools","info"] + sketch_paths
        # cmd = " ".join(cmdlist)
        return cmd
    

    def sketch_check(self, path=None) -> bool:
        '''Check that kmc databases associated with full path exist and are not empty'''
        # If path is unknown take it from the sfp
        if not path:
            path = self.sfp.full
        
        # First check if files exist and are non-empty
        if not (os.path.exists(path + ".kmc_pre") and 
                os.path.exists(path + ".kmc_suf") and 
                os.stat(path + ".kmc_suf").st_size != 0 and 
                os.stat(path + ".kmc_pre").st_size != 0):
            return False
        
        # Then check if we have a valid cardinality stored
        if path in self.speciesinfo.cardkey:
            return self.speciesinfo.cardkey[path] > 0
        
        return True  # Files exist but no cardinality stored yet
        
    def remove_sketch(self, delete_me:str=None):
        '''Delete kmc database files for the associated sketch object'''
        if not delete_me:
            sketchname = self.sfp.full
        else:
            sketchname = delete_me
        if self.kval == 0:
            sketchname = self.sfp.full.replace("{}","*")
        try:
            if '*' in sketchname:
                for f in glob.glob(sketchname):
                    for ext in [".kmc_pre", ".kmc_suf"]:
                        if os.path.exists(f + ext):
                            os.remove(f + ext)
            else:
                for ext in [".kmc_pre", ".kmc_suf"]:
                    if os.path.exists(sketchname + ext):
                        os.remove(sketchname + ext)
        except FileNotFoundError:
            pass

    def _leaf_command(self,tmpdir) -> str:
        '''Command string to produce the kmc database files from a fasta based on the information used to initiate the sketch obj'''
        kval_str=str(self.kval)
        if self.kval == 0:
            kval_str="{}"
        tmpkdir=os.path.join(tmpdir,"k"+kval_str)
        os.makedirs(tmpkdir, exist_ok=True)
        tstring = ''
        if self.experiment["nthreads"] > 0:
            tstring = ' -t'+ str(self.experiment['nthreads'])
        cmdlist = ['kmc -hp' + tstring,
        ' -ci1 -cs2','-k' + kval_str,
        canon_command(canon=self.experiment['canonicalize'], tool="kmc"),
        '-fm', self.sfp.ffiles[0], self.sfp.full, tmpkdir]
        cmd = " ".join(cmdlist)
        return cmd

    def _union_command(self) -> str:
        '''Command string to produce unions of the kmc database files based on the information used to initiate the sketch obj'''
        complex_input = "INPUT: \n"
        inputn=1
        tstring = ''
        if self.experiment["nthreads"] > 0:
            tstring = ' -t'+ str(self.experiment['nthreads'])

        for sketch in self._presketches:
            complex_input = complex_input + f"input{inputn} = {sketch} -ci1   \n"
            inputn+=1
        complex_input = complex_input + f"OUTPUT:\n{self.sfp.full} = " + " + ".join([f"input{i+1}" for i in range(inputn-1)])
        cmdlist = [f'echo -e "{complex_input}"',"|","kmc_tools","-hp ",tstring, "complex", "/dev/stdin"]
        cmd = " ".join(cmdlist)
        return cmd
