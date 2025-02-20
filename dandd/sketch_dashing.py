from __future__ import annotations
import os
import subprocess
import csv
import glob
from dandd.sketch_base import SketchObj
from dandd.utils import canon_command

DASHINGLOC="/home/jbonnie1/lib/dashing_working/dashing"
# DASHINGLOC = "dashing"
   
class DashSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches)
    
    def card_command(self, sketch_paths=[]) -> str:
        '''Return the cardinality command for the provided dashing sketch paths'''
        if len(sketch_paths) == 0:
            if self.kval == 0:
                return
            sketch_paths=[self.sketch]
        cmdlist = [DASHINGLOC,"card", "--presketched"]+ sketch_paths #+ verbose
        #, "-p10"
        # , f"-p{self.experiment['nthreads']}"
        cmd = " ".join(cmdlist)
        return cmd

    def parse_card(self, proc):
        '''Parse the cardinality streaming from standard out for the dashing card command'''
        for card in csv.DictReader(proc.stdout.splitlines(),delimiter='\t'):
            self.speciesinfo.cardkey[card['#Path']] = float(card['Size (est.)'])

    def sketch_check(self, path=None) -> bool:
        '''Check that dashing sketch at full path exists and is not empty'''
        if not path:
            path=self.sfp.full
        # if self.experiment["lowmem"]:
        #     # print(self.check_cardinality())
        #     if self.check_cardinality() > 0:
        #         return True
        if os.path.exists(path) and os.stat(path).st_size != 0:
            return True
        else:
            return False
    
    def remove_sketch(self, delete_me:str=None):
        '''Delete intermediate sketches in batches'''
        if not delete_me:
            sketchname= self.sfp.full
        else:
            sketchname = delete_me
        if self.kval == 0:
            sketchname = self.sfp.full.replace("{}","*")
        try:
            for f in glob.glob(sketchname):
                os.remove(f)
        except FileNotFoundError:
            # print(f"{sketchname} has already been removed.")
            pass
    
    def _leaf_command(self, tmpdir) -> str:
        '''Command string to produce the sketch from a fasta based on the information used to initiate the dashing sketch obj'''
        str_kval=str(self.kval)
        file_source= [self.sfp.ffiles[0]]
        if self.kval == 0:
            str_kval="{}"
        
        cmdlist = [DASHINGLOC, "sketch", 
        canon_command(self.experiment['canonicalize'], "dashing"),
        "-k" + str_kval, 
        "-S",str(self.experiment['registers']),
        #  f"-p{self.experiment['nthreads']}",
         "--prefix", str(self.sfp.dir)]
        
        cmd = " ".join(cmdlist+file_source)
        return cmd
    
    def _union_command(self) -> str:
        '''Returns bash command to create a union sketch'''
        cmdlist = [DASHINGLOC, "union", #f"-p{self.experiment['nthreads']}",
        "-z -o", str(self.sfp.full)] + self._presketches #+ verbose
        cmd = " ".join(cmdlist)
        return cmd
        
