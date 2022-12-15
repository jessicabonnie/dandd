from species_specifics import SpeciesSpecifics
import os
import hashlib
import warnings
import subprocess
import csv
import sys
import tempfile
import shutil

DASHINGLOC="dashing" #"/scratch16/blangme2/jessica/lib/dashing/dashing"

def blake2b(fname):
    hash_blake2b = hashlib.blake2b()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_blake2b.update(chunk)
    return hash_blake2b.hexdigest()

def canon_command(canon:bool, tool='dashing'):
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
        # fasta_list="\n".join(self._sfp.ffiles)
    if experiment["tool"] == "dashing":
        progeny_dir=os.path.join(sketchdir, "ngen" + str(1),"k"+ str(kval))
        cmdlist = [ DASHINGLOC, "sketch", 
    canon_command(experiment['canonicalize'], "dashing"),
    "-k" + str(kval), 
    "-S",str(experiment['registers']),
        "-p10","--prefix", progeny_dir, "--paths", "/dev/stdin"]
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
        #self.baseold = self.nameSketch(speciesinfo=speciesinfo, kval=kval)
        self.base =self.assign_base(speciesinfo=speciesinfo, kval=kval, registers=experiment['registers'], canonicalize=experiment['canonicalize'])
        self.dir = os.path.join(speciesinfo.sketchdir, "ngen" + str(self.ngen),"k"+ str(kval))
        self.full = os.path.join(self.dir, self.base)+ self._get_ext(experiment['tool'])
        #self.registers = speciesinfo.registers
        os.makedirs(self.dir, exist_ok=True) 
    #TODO: WHY IS ngen10 happening without the sketch
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
        
    def hashsum(self, speciesinfo:SpeciesSpecifics):
        if self.ngen == 1:
            output= blake2b(self.ffiles[0])
        else:
            sum = int("0",16)
            for fasta in self.files:
                sum += int(speciesinfo.fastahex[fasta],16)
            output= hex(sum)
        return output
        
    def assign_base(self, speciesinfo:SpeciesSpecifics, kval:int, registers:int, canonicalize:bool):
        fnames_key=''.join(self.files)
        if fnames_key not in speciesinfo.fastahex.keys():
            speciesinfo.fastahex[fnames_key]= self.hashsum(speciesinfo)   
            stored_val=speciesinfo.fastahex[fnames_key]
        else:
            checkval=self.hashsum(speciesinfo)
            stored_val=speciesinfo.fastahex[fnames_key]
            if checkval != stored_val:
               raise RuntimeError(f"Checksum does not match stored value for {fnames_key}: {checkval}, {stored_val}")
        if self.ngen == 1:
            sketchbase=self.files[0] + ".w." + str(kval) + ".spacing." + str(registers)
        else:
            suffix = str(registers) + "n" + str(self.ngen) + "k" + str(kval)
            if not canonicalize:
                suffix=suffix+'nc'
            sketchbase = stored_val[:15] + "_" + suffix
        # store information relating to this basename to be given to user later as table or obj
        info = [self.files,self.ngen, kval, registers ]
        if sketchbase not in speciesinfo.sketchinfo.keys():
            speciesinfo.sketchinfo[sketchbase] = info
        else:
            stored_info=speciesinfo.sketchinfo[sketchbase]
            for index in range(len(stored_info)):
              if stored_info[index] != info[index]:
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
    
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        self.kval = kval
        self.experiment=experiment
        self.sketch = None
        self.cmd = None
        self._sfp = sfp
        self._presketches = presketches
        self.create_sketch()
        self.card = self.check_cardinality(speciesinfo=speciesinfo)
        self.delta_pos = self.card/self.kval
    
    def __lt__(self, other):
        # lt = less than
        return self.delta_pos < other.delta_pos
    def __gt__(self, other):
        # gt = greater than
        return self.delta_pos > other.delta_pos
    def __repr__(self):
        return f"['sketch loc: {self.sketch}', k: {self.kval}, pos delta: {self.delta_pos}, cardinality: {self.card}, command: {self.cmd}  ]"
    
    def sketch_usable(self) -> bool:
        raise NotImplementedError
    
    def sketch_check(self)->bool:
        raise NotImplementedError("Subclass needs to define this.")
        
    def leaf_command(self, tmpdir) -> str:
        raise NotImplementedError("Subclass needs to define this.")
    def command_check(self):
        pass

    def leaf_sketch(self, just_do_it=False):
        ''' If leaf sketch file exists, record the command that would have been used. If not run the command and store it.'''
        tmpdir=tempfile.mkdtemp()
        cmd = self.leaf_command(tmpdir=tmpdir)
        print(f"Just do it: {just_do_it}")
        if self.experiment['debug']:
            print(cmd)
        if just_do_it:
            print("Just Do It.")
            subprocess.call(cmd, shell=True)
            print("Just did it.")
            self.cmd=cmd
        elif not self.sketch_check():
            print("The sketch file {0} either doesn't exist or is empty".format(self._sfp.full))
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
            
            ##TODO pass back the command so that it can be done in parallel?
        else:
            self.cmd = cmd

        shutil.rmtree(tmpdir)
        #print(self.cmd)

    ##TODO: need to separate process of identifying and saving hashkey so a badly formed sketch/db doesn't get saved in the key... or we catch the associated error and delete the sketch file and try again with a new one 
    #     filename =os.path.basename(filename)
    
    def union_command(self)->str:
        raise NotImplementedError("Subclass needs to define this.")

    
    def union_sketch(self, just_do_it=False):
        ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmd = self.union_command()
        if just_do_it:
            subprocess.call(cmd, shell=True)
            self.cmd = cmd
            # subprocess.call(cmd, shell=True)
        elif not self.sketch_check():
            print(f"The sketch file {self._sfp.full} either doesn't exist or is empty")
            #self.cmd=subprocess.run(cmdlist)
            subprocess.call(cmd, shell=True)
            # subprocess.call(cmd, shell=True)
            self.cmd = cmd
        
        self.cmd = cmd
            #" ".join(cmdlist)
        if self.experiment['debug']:
            print(self.cmd)
    ##TODO: make leaf sketch and union sketch private
        
    def create_sketch(self, just_do_it=False):
        ''' If sketch file exists, assign path to self.sketch and return path. 
            If not create sketch, assign, and then return path.'''
        if self._sfp.ngen == 1:
            print("INSIDE")
            self.leaf_sketch(just_do_it=just_do_it)
        elif self._sfp.ngen > 1:
            self.union_sketch(just_do_it=just_do_it)
        else:
            raise RuntimeError("For some reason you are trying to sketch an empty list of files. Don't do that.")
        
        self.sketch = self._sfp.full
        return self.sketch
    

    def individual_card(self, speciesinfo:SpeciesSpecifics):
        '''Run cardinality for an individual sketch or database. Add it to a dictionary {path:value}'''
        raise NotImplementedError("Subclass needs to define this.")
        
    def check_cardinality(self, speciesinfo: SpeciesSpecifics):
        print("inside check_cardinality")
        if self._sfp.full not in speciesinfo.cardkey.keys() or speciesinfo.cardkey[self._sfp.full] == 0:
            self.individual_card(speciesinfo=speciesinfo)   
        return float(speciesinfo.cardkey[self.sketch])

class DashSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches)
    
    def card_command(self, sketch_paths) -> str:
        cmdlist = [DASHINGLOC,"card", "--presketched", "-p10"]+ sketch_paths
        cmd = " ".join(cmdlist)
        return cmd
    
    
    def view_test_command(self) -> str:
        '''Check to make sure that the supposedly existing sketch is viewable / not malformed. NOT CURRENTLY IN USE'''
        cmdlist = [DASHINGLOC,"view", self._sfp.full]
        cmd = " ".join(cmdlist)
        return cmd

    def individual_card(self, speciesinfo: SpeciesSpecifics):
        '''Run cardinality for an individual sketch. Add it to a dictionary {path:value}'''
        cmd = self.card_command([self._sfp.full])
        
        if self.experiment['debug']:
            print(cmd)
        #print(cmd)
        try:
            print("trying card")
            proc=subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, check=True,universal_newlines=True)
            print("tried")
            # for card in csv.DictReader(card_lines, delimiter='\t'):
            #     speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
            # readme=False
            # print(card_lines)
            # for card in card_lines:
            #     llist=card.split()
            #     print(llist)
            #     if readme:
            #         llist=card.split()
            #         print(llist)
            #         speciesinfo.cardkey[llist[0]] = llist[1]
            #         readme=True

        except subprocess.CalledProcessError:
            print("failed")
            # warnings.warn(message=f"{self._sfp.full} cannot be created. Will attempt to remove and recreate component sketches.", category=RuntimeWarning)
            self.create_sketch(just_do_it=True)
            print("recreated sketch")
            proc=subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, universal_newlines=True)
            print("ran card again")
            # card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout
        finally:
            print("I think I have a working sketch")
            for card in csv.DictReader(proc.stdout.splitlines(),delimiter='\t'):
                speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
        #return super().individual_card(speciesinfo, debug)

    def parse_card(self,speciesinfo, proc):
        for card in csv.DictReader(proc.stdout.splitlines(),delimiter='\t'):
                speciesinfo.cardkey[card['#Path']] = card['Size (est.)']

    def sketch_check(self) -> bool:
        '''Check that sketch at full path exists and is not empty'''
        
        if os.path.exists(self._sfp.full) and os.stat(self._sfp.full).st_size != 0:
            return True
        else:
            return False
            # and os.stat(self._sfp.full).st_size != 0:
        #     if self.sketch_usable():
        #         return True
        # else:
        #     return False
    
    def sketch_usable(self)-> bool:
        '''Not currently in use'''
        try:
            code=subprocess.run(self.view_test_command(), shell=True,check=True, stdout=subprocess.DEVNULL,start_new_session=True).returncode
            if code != 0:
                print("INSIDE")
                warnings.warn(message=f"{self._sfp.full} cannot be created. Will attempt to remove and recreate component sketches.", category=UserWarning)
                return False
        except:
            warnings.warn(message=f"{self._sfp.full} cannot be created. Will attempt to remove and recreate component sketches.", category=UserWarning)
            return False
        else:
            return True
    
    def leaf_command(self, tmpdir) -> str:
        '''Command string to produce the sketch from a fasta based on the information used to initiate the sketch obj'''
        cmdlist = [DASHINGLOC, "sketch", 
        canon_command(self.experiment['canonicalize'], "dashing"),
        "-k" + str(self.kval), 
        "-S",str(self.experiment['registers']),
         "-p10","--prefix", str(self._sfp.dir),
          self._sfp.ffiles[0]]
                #    os.path.join(speciesinfo.inputdir, sfp.files[0])]
        cmd = " ".join(cmdlist)
        return cmd

    
    def union_command(self) -> str:
        '''Returns bash command to create a union sketch'''
        cmdlist = [DASHINGLOC, "union", "-p 10 ","-z -o", str(self._sfp.full)] + self._presketches
        cmd = " ".join(cmdlist)
        return cmd
        
    # def union_sketch(self):
    #     ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
    #     cmdlist = [DASHINGLOC, "union", "-p 10 ","-z -o", str(self._sfp.full)] + self._presketches
    #     cmd = " ".join(cmdlist)
    #     #print("SIZE of {0} is {1}".format(self._sfp.full, os.stat(sfp.full)))
    #     if (not os.path.exists(self._sfp.full)) or os.stat(self._sfp.full).st_size == 0:
    #         # print("The sketch file {0} either doesn't exist or is empty".format(sfp.full))
    #         #self.cmd=subprocess.run(cmdlist)
    #         subprocess.call(cmd, shell=True)
    #         self.cmd=cmd
    #     else:
    #         self.cmd = cmd
    #         #" ".join(cmdlist)
    #     if self.experiment['debug']:
    #         print(self.cmd)

class KMCSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[]):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches)

    def individual_card(self, speciesinfo: SpeciesSpecifics):
        cmdlist = ["kmc_tools","info",self.sketch] #this should be self._sfp.full
        cmd = " ".join(cmdlist)
        if self.experiment['debug']:
            print(cmd)
        card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
        for line in card_lines:
            print(line)
            key, value = line.strip().split(':')
            print (f"{key},{value}")
            if key.strip() == 'total k-mers':
                speciesinfo.cardkey[self.sketch] = value.strip()
        #return super().individual_card(speciesinfo, debug)
    
    def sketch_check(self) -> bool:
        if (os.path.exists(self._sfp.full +".kmc_pre")) and (os.path.exists(self._sfp.full +".kmc_suf")) and os.stat(self._sfp.full +".kmc_suf").st_size != 0:
            return True
        else:
            return False

    def leaf_command(self,tmpdir) -> str:
        cmdlist = ['kmc -t'+ str(self.experiment['nthreads']),
        '-ci1 -cs2',f'-k{str(self.kval)}',
        canon_command(canon=self.experiment['canonicalize'], tool="kmc"),
        '-fm', self._sfp.ffiles[0], self._sfp.full, tmpdir]
        cmd = " ".join(cmdlist)
        return cmd

    def union_command(self) -> str:
        complex_input = "INPUT: \n"
        inputn=1
        for sketch in self._presketches:
            complex_input = complex_input + f"input{inputn} = {sketch} -ci1   \n"
            inputn+=1
        complex_input = complex_input + f"OUTPUT:\n{self._sfp.full} = " + " + ".join([f"input{i+1}" for i in range(inputn-1)])
        cmdlist = [f'echo -e "{complex_input}"',"|","kmc_tools","-t"+ str(self.experiment['nthreads']), "complex", "/dev/stdin"]
        cmd = " ".join(cmdlist)
        return cmd
