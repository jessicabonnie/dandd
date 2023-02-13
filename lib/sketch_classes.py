from species_specifics import SpeciesSpecifics
import os
import hashlib
import subprocess
import csv
import tempfile
import shutil

# This assumes that the command for dashing has been aliased
DASHINGLOC="dashing" 

def blake2b(fname):
    '''Create a blake2b hexsum from a file'''
    hash_blake2b = hashlib.blake2b()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_blake2b.update(chunk)
    return hash_blake2b.hexdigest()

def canon_command(canon:bool, tool='dashing'):
    '''Determine what should be added to sketching command when not canonicalizing
    '''
    outstr=''
    if not canon:
        if tool == 'dashing':
            outstr='--no-canon'
        if tool == 'kmc':
            outstr='-b'
    return outstr


def parallel_progeny_command(sketchdir, kval:int, experiment) -> str:
    '''Produce the necessary leaf sketches for a node further up the tree in a single batch command to dashing/kmc
    NOTE: NOT CURRENTLY USED'''
        # fasta_list="\n".join(self.sfp.ffiles)
    if experiment["tool"] == "dashing":
        progeny_dir=os.path.join(sketchdir, "ngen" + str(1),"k"+ str(kval))
        cmdlist = [ DASHINGLOC, "sketch", 
    canon_command(experiment['canonicalize'], "dashing"),
    "-k" + str(kval), 
    "-S",str(experiment['registers']),
        #f"-p{experiment['nthreads']}",
        "--prefix", progeny_dir, "--paths", "/dev/stdin"]
            #    os.path.join(speciesinfo.inputdir, sfp.files[0])]
    else:
        cmdlist=[]
        raise NotImplementedError("tools other than dashing not yet implemented with parallel argument")
    cmd = " ".join(cmdlist)
    
    return cmd

class SketchFilePath:
    '''An object to prepare sketch and union naming and directory location
    filenames: list of filenames that will be used to make the sketch
    kval: kvalue for the sketch
    speciesinfo: SpeciesSpecifics object for the species
    prefix: currently unused tag for filenames to differentiate between runs
    '''
    def __init__(self, filenames: list, kval: int, speciesinfo: SpeciesSpecifics, experiment:dict, prefix=None ):
        self.ffiles = filenames
        self.files= [os.path.basename(f) for f in self.ffiles]
        self.files.sort()
        self.ngen = len(filenames)
        self.dir = os.path.join(speciesinfo.sketchdir, "ngen" + str(self.ngen),"k"+ str(kval))
        #self.baseold = self.nameSketch(speciesinfo=speciesinfo, kval=kval)
        self.base =self._assign_base(speciesinfo=speciesinfo, kval=kval, registers=experiment['registers'], canonicalize=experiment['canonicalize'], tool=experiment['tool'], safety=experiment['safety'])
        self.relative = os.path.join("ngen" + str(self.ngen),"k"+ str(kval),self.base)+ self._get_ext(experiment['tool'])
        # self.full = os.path.join(speciesinfo.sketchdir, self.relative)
        self.full = os.path.join(self.dir, self.base)+ self._get_ext(experiment['tool'])
        #self.registers = speciesinfo.registers
        os.makedirs(self.dir, exist_ok=True) 
    
    #WHAT K/D are reported when tree node is reported.
    def __repr__(self):
        return f"{self.__class__.__name__}[basename: {self.base}, 'fullpath inputs: {self.ffiles}', ngen: {self.ngen}, dir: {self.dir}, fullpath: {self.full} ]"
    def _get_ext(self,tool)->str:
        if tool == 'dashing':
            ext='.hll'
        elif tool == 'kmc':
            ext=''
        else:
            raise ValueError("is there another option for tool other than kmc or dashing?")
        return ext
        
    def _hashsum(self, speciesinfo:SpeciesSpecifics):
        '''Calculate the blake2b hexsum of an individual fastas or sum the hexsums of component fastas to create hexidecimal identifiers for combinations of fastas'''
        if self.ngen == 1:
            output= blake2b(self.ffiles[0])
        else:
            sum = int("0",16)
            for fasta in self.files:
                sum += int(speciesinfo.fastahex[fasta],16)
            output= hex(sum)
        return output
        
    def _assign_base(self, speciesinfo:SpeciesSpecifics, kval:int, registers:int, canonicalize:bool, tool:str, safety=False) -> str:
        '''determine the base file name for the sketch using the properties that will be used to generate it'''
        fnames_key=''.join(self.files)
        # if the key (made by joining the ingredient filenames) isn't already in the fastahex dictionary mapping the combination of those files to a hexsum, calculate that hexsum and add it to the fastahex key
        if fnames_key not in speciesinfo.fastahex.keys():
            speciesinfo.fastahex[fnames_key]= self._hashsum(speciesinfo)   
            stored_val=speciesinfo.fastahex[fnames_key]
        # if the key is there, calculate what we expect the hashsum value to be based on the hexes of the components -- this is just to confirm that nothing has gotten confused somehow
        else:
            stored_val=speciesinfo.fastahex[fnames_key]
            if safety:
                checkval=self._hashsum(speciesinfo)
                if checkval != stored_val:
                    raise RuntimeError(f"Checksum does not match stored value for {fnames_key}: {checkval}, {stored_val}")
        # if the sketch is of a single input file, dashing will insist on naming it something specific, so we will use that base as a name for both dashing and kmc to make life easier
        if self.ngen == 1:
            if tool == 'dashing':
                sketchbase=self.files[0] + ".w." + str(kval) + ".spacing." + str(registers)
            else:
                sketchbase=self.files[0] + "_k" + str(kval)
                if not canonicalize:
                    sketchbase=sketchbase+'nc'
        # if sketch is of a combination of sketches, add a tag that will differentiate it from other combinations of the same sketches composed using different register counts or kvalues. Also guard against the low probability chance that there are overlapping hexsums of different numbers of fasta inputs
        else:
            suffix = str(registers) + "n" + str(self.ngen) + "k" + str(kval)
            if not canonicalize:
                suffix=suffix+'nc'
            sketchbase = stored_val[:15] + "_" + suffix
        # store information relating to this basename to be given to user later as table or obj
        info = {"sketchbase": sketchbase, "files": self.files, "ngen": self.ngen, "kval": kval, "registers": registers }
        if sketchbase not in speciesinfo.sketchinfo.keys():
            speciesinfo.sketchinfo[sketchbase] = info
        elif safety:
            stored_info=speciesinfo.sketchinfo[sketchbase]
            # check to make sure all sketchinfo values match what is stored
            for key in stored_info.keys():
              if stored_info[key] != info[key]:
                  raise RuntimeError(f"Duplicate keys but not duplicate values: {sketchbase}: (1) {stored_info}, (2) {info}")
        return sketchbase


class SketchObj(object):
    ''' A sketchobject.
        Inputs:
            kval: kmer length
            sfp: SketchFilePath object
            speciesinfo: SpeciesSpecifics object containing path information for the species/tag
            presketches: list of sketches to use in union sketching

        Properties:
            kval = the kvalue used to construct the sketch
            sketch = the location of the sketch file
            cmd = the command used to create the sketch
            card = the cardinality of the sketch
            delta_pos = possible delta value (to be compared to other ks of the same group of files)
    '''
    
    def __init__(self, kval: int, sfp: SketchFilePath, speciesinfo: SpeciesSpecifics, experiment: dict, presketches=[]):
        self.kval = kval
        self.sketch = None
        self.cmd = None
        self.sfp = sfp
        self.speciesinfo=speciesinfo
        experiment['baseset'].add(sfp.base)
        self.experiment=experiment
        self._presketches = presketches
        self.create_sketch()
        self.card = self.check_cardinality()
        self.delta_pos = self.card/self.kval
    
    def __lt__(self, other):
        # lt = less than
        return self.delta_pos < other.delta_pos
    def __gt__(self, other):
        # gt = greater than
        return self.delta_pos > other.delta_pos
    def __repr__(self):
        return f"['sketch loc: {self.sketch}', k: {self.kval}, pos delta: {self.delta_pos}, cardinality: {self.card}, command: {self.cmd}  ]"
    
    def sketch_check(self)->bool:
        raise NotImplementedError("Subclass needs to define this.")
        
    def _leaf_command(self, tmpdir) -> str:
        raise NotImplementedError("Subclass needs to define this.")
    def command_check(self):
        pass

    def _leaf_sketch(self, just_do_it=False):
        ''' If leaf sketch file exists, record the command that would have been used. If not run the command and store it.'''
        tmpdir=tempfile.mkdtemp()
        cmd = self._leaf_command(tmpdir=tmpdir)
        if self.experiment['debug']:
            print(cmd)
        if just_do_it:
            if self.experiment['verbose']:
                print("Due to issues with leaf sketch/db file, we will Just Do It. (It=Sketch or Build Again)")
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
        elif not self.sketch_check():
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
            
            ##TODO pass back the command so that it can be done in parallel?
        else:
            self.cmd = cmd

        shutil.rmtree(tmpdir)
        #print(self.cmd)
    #     filename =os.path.basename(filename)
    
    def _union_command(self)->str:
        raise NotImplementedError("Subclass needs to define this.")

    
    def _union_sketch(self, just_do_it=False):
        ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmd = self._union_command()
        if just_do_it:
            subprocess.call(cmd, shell=True)
            self.cmd = cmd
        elif not self.sketch_check():
            # print(f"The sketch file {self.sfp.full} either doesn't exist or is empty")
            subprocess.call(cmd, shell=True)
            self.cmd = cmd
        
        self.cmd = cmd
        if self.experiment['debug']:
            print(self.cmd)
        
    def create_sketch(self, just_do_it=False):
        ''' If sketch file exists, assign path to self.sketch and return path. 
            If not create sketch, assign, and then return path.'''
        if self.sfp.ngen == 1:
            self._leaf_sketch(just_do_it=just_do_it)
        elif self.sfp.ngen > 1:
            self._union_sketch(just_do_it=just_do_it)
        else:
            raise RuntimeError("For some reason you are trying to sketch an empty list of files. Don't do that.")
        
        self.sketch = self.sfp.full
        return self.sketch
    
    def card_command(self, sketch_paths=[]):
        raise NotImplementedError("Subclass must define this.")
    def parse_card(self, proc):
        raise NotImplementedError("Subclass must define this.")

    def individual_card(self):
        '''Run cardinality for an individual sketch or database. Add it to a dictionary {path:value}'''
        cmd = self.card_command()#[self.sfp.full])
        if self.experiment['debug']:
            print(cmd)
        try:
            proc=subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, check=True,universal_newlines=True)
        except subprocess.CalledProcessError:
            # warnings.warn(message=f"{self.sfp.full} cannot be created. Will attempt to remove and recreate component sketches.", category=RuntimeWarning)
            self.create_sketch(just_do_it=True)
            # print("recreated sketch")
            proc=subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, universal_newlines=True)
        finally:
            self.parse_card(proc=proc)
        
    def check_cardinality(self, delay=False):
        '''Check whether the cardinality of sketch/db is stored in the cardkey, if not run a card command for the sketch and store it. NOTE: delay argument not implemented relates to idea of a way to batch the card command by attatching the sketch to a card0 sketch list attached to the SpeciesSpecifics object'''
        
        if self.sfp.full not in self.speciesinfo.cardkey.keys():
            if self.experiment['verbose']:
                print(f"Sketch File Path not in cardkey: {self.sfp.full}")
        if self.sfp.full not in self.speciesinfo.cardkey.keys() or self.speciesinfo.cardkey[self.sfp.full] == 0:
            self.individual_card()
        return float(self.speciesinfo.cardkey[self.sketch])
    
    # def summarize(self):
    #     outdict = {"kval":self.kval, "card":self.card, "delta_pos":self.delta_pos}
    #     return outdict
class DashSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches)
    
    def card_command(self, sketch_paths=[]) -> str:
        if len(sketch_paths) == 0:
            sketch_paths=[self.sketch]
        cmdlist = [DASHINGLOC,"card", "--presketched"]+ sketch_paths
        #, "-p10"
        # , f"-p{self.experiment['nthreads']}"
        cmd = " ".join(cmdlist)
        return cmd

    def parse_card(self, proc):
        '''Parse the cardinality streaming from standard out for the card command'''
        for card in csv.DictReader(proc.stdout.splitlines(),delimiter='\t'):
                self.speciesinfo.cardkey[card['#Path']] = card['Size (est.)']

    def sketch_check(self) -> bool:
        '''Check that sketch at full path exists and is not empty'''
        
        if os.path.exists(self.sfp.full) and os.stat(self.sfp.full).st_size != 0:
            return True
        else:
            return False
    
    
    def _leaf_command(self, tmpdir) -> str:
        '''Command string to produce the sketch from a fasta based on the information used to initiate the sketch obj'''
        cmdlist = [DASHINGLOC, "sketch", 
        canon_command(self.experiment['canonicalize'], "dashing"),
        "-k" + str(self.kval), 
        "-S",str(self.experiment['registers']),
        #  f"-p{self.experiment['nthreads']}",
         "--prefix", str(self.sfp.dir),
          self.sfp.ffiles[0]]
                #    os.path.join(speciesinfo.inputdir, sfp.files[0])]
        cmd = " ".join(cmdlist)
        return cmd
    
    def _union_command(self) -> str:
        '''Returns bash command to create a union sketch'''
        cmdlist = [DASHINGLOC, "union", #f"-p{self.experiment['nthreads']}",
        "-z -o", str(self.sfp.full)] + self._presketches
        cmd = " ".join(cmdlist)
        return cmd
        


class KMCSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches)

    
    def parse_card(self, proc):
        for line in proc.stdout.splitlines():
            key, value = line.strip().split(':')
            # print (f"{key},{value}")
            if key.strip() == 'total k-mers':
                self.speciesinfo.cardkey[self.sketch] = value.strip()

    def card_command(self, sketch_paths:list=[]) -> str:
        if len(sketch_paths) == 0:
            sketch_paths=[self.sketch]
        cmdlist = ["kmc_tools","info"] + sketch_paths
        cmd = " ".join(cmdlist)
        return cmd
    def sketch_check(self) -> bool:
        if (os.path.exists(self.sfp.full +".kmc_pre")) and (os.path.exists(self.sfp.full +".kmc_suf")) and os.stat(self.sfp.full +".kmc_suf").st_size != 0:
            return True
        else:
            return False

    def _leaf_command(self,tmpdir) -> str:
        cmdlist = ['kmc -t'+ str(self.experiment['nthreads']),
        '-ci1 -cs2',f'-k{str(self.kval)}',
        canon_command(canon=self.experiment['canonicalize'], tool="kmc"),
        '-fm', self.sfp.ffiles[0], self.sfp.full, tmpdir]
        cmd = " ".join(cmdlist)
        return cmd

    def _union_command(self) -> str:
        complex_input = "INPUT: \n"
        inputn=1
        for sketch in self._presketches:
            complex_input = complex_input + f"input{inputn} = {sketch} -ci1   \n"
            inputn+=1
        complex_input = complex_input + f"OUTPUT:\n{self.sfp.full} = " + " + ".join([f"input{i+1}" for i in range(inputn-1)])
        cmdlist = [f'echo -e "{complex_input}"',"|","kmc_tools","-t"+ str(self.experiment['nthreads']), "complex", "/dev/stdin"]
        cmd = " ".join(cmdlist)
        return cmd
