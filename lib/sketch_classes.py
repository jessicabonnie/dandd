from species_specifics import SpeciesSpecifics
import os
import hashlib
import warnings
import subprocess
import csv

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
        if tool == 'kmc':
            ext=''
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


    # def assign_hash_string(self, filename: str, speciesinfo: SpeciesSpecifics, length: int):
    #     '''If the basename of filename already exists in hash/dict, look it up; otherwise create one and add to hash/dict using the full path as the key'''
    #     ##TODO: need to separate process of identifying and saving hashkey so a badly formed sketch/db doesn't get saved in the key... or we catch the associated error and delete the sketch file and try again with a new one 
    #     filename =os.path.basename(filename)
    #     if filename in speciesinfo.hashkey.keys():
    #         return speciesinfo.hashkey[filename]
    #     else:
    #         alphanum=hashlib.md5(filename.encode()).hexdigest()
    #         trunc=alphanum[:length]
    #         if trunc in speciesinfo.hashkey.values():
    #             #warnings.warn("Hashvalue " + trunc + " has 2 keys!! " + filename + " will be assigned to a longer hash.")
    #             return self.assign_hash_string(filename, speciesinfo, length=length+1)
    #         else:
    #             speciesinfo.hashkey[filename] = trunc
    #             return trunc
    
    ##TODO: implement naming function for kmc
    # def nameSketch(self, speciesinfo: SpeciesSpecifics, kval: int):
    #     '''Determine what the name of a sketch is/will be'''
    #     if self.ngen > 1:
            
    #         filehashes = [self.assign_hash_string(filename=onefile, speciesinfo=speciesinfo, length=1) for onefile in self.files]
    #         filehashes.sort()
    #         outfile_prefix = speciesinfo.tag + "_" + "_".join(filehashes) + "_k" + str(kval) + "_r" + str(speciesinfo.registers)
    #         if len(outfile_prefix) >30:
    #             alphanum=hashlib.md5(outfile_prefix.encode()).hexdigest()
    #             outfile_prefix = alphanum
    #         outfile=outfile_prefix + ".hll"
    #     else:
    #         outfile=os.path.basename(self.files[0]) + ".w." + str(kval) + ".spacing." + str(speciesinfo.registers) + ".hll"
    #     return outfile    

##TODO create superclasses for precise vs estimate classes which have same function names
class SketchObj:
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
    
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[], delay=False, debug=True):
        self.kval = kval
        self.experiment=experiment
        self.sketch = None
        self.cmd = None
        self._sfp = sfp
        self._presketches = presketches
        self.create_sketch(debug=debug)
        self.card = self.check_cardinality(speciesinfo=speciesinfo, debug=debug)
        self.delta_pos = self.card/self.kval
    
    def __lt__(self, other):
        # lt = less than
        return self.delta_pos < other.delta_pos
    def __gt__(self, other):
        # gt = greater than
        return self.delta_pos > other.delta_pos
    def __repr__(self):
        return f"['sketch loc: {self.sketch}', k: {self.kval}, pos delta: {self.delta_pos}, cardinality: {self.card}, command: {self.cmd}  ]"
    
    '''TODO: implement kmc
    CODE from benchmarking to create individual databases:
    if [ $approach == 'kmc' ]; then
    outsketch=${apout}/${fasta}.kmc
    /usr/bin/time -o ${outprefix}.out -v sh -c "kmc -v -t${nthreads} -k${kval} -ci1 -fm ${datadir}/${fasta} ${apout}/${fasta}.kmc ${outdir}/kmc > ${cardloc}"
    card=$(grep "No. of unique counted k-mers" ${cardloc} | awk '{print $NF}')
    '''
    
    def sketch_check(self):
        if (os.path.exists(self._sfp.full +".kmc_pre")) and (os.path.exists(self._sfp.full +".kmc_suf")) and os.stat(self._sfp.full +".kmc_suf").st_size == 0:
            return True
        else:
            return False
        
    def leaf_command(self):
        pass
    
    def leaf_sketch(self, debug=False, delay=False, nthreads=10):
        ''' If leaf sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmd = self.leaf_command()
        if debug:
                print(cmd)
        if not self.sketch_check():
            print("The sketch file {0} either doesn't exist or is empty".format(self._sfp.full))
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
            ##TODO Check if call raises error / returns 0
            ##TODO pass back the command so that it can be done in parallel?
        else:
            self.cmd = cmd
        #print(self.cmd)

    ##TODO: need to separate process of identifying and saving hashkey so a badly formed sketch/db doesn't get saved in the key... or we catch the associated error and delete the sketch file and try again with a new one 
    #     filename =os.path.basename(filename)
    

    '''TODO: implement kmc
    CODE from benchmarking to create union from database list
 
  #Create an instruction file for kmc_tools complex
  echo "INPUT:" > ${apout}/complex_union.txt
  echo ${sketched[@]} | sed 's/ /\n/g' | awk '{print "input"NR" = ",$0, "-ci1"}'>> ${apout}/complex_union.txt
  echo "OUTPUT:" >> ${apout}/complex_union.txt
  string="${fullunion} = input1"
  for i in $(seq 2 $nfasta); do string="${string} + input$i"; done
  echo ${string} >> ${apout}/complex_union.txt
  
  #Benchmarch the process
  # TODO : update to use info command 

  cmd="kmc_tools -t${nthreads} complex ${apout}/complex_union.txt" 
  echo ${cmd}
  /usr/bin/time -o ${timeout} -v sh -c "${cmd}"
  kmc_tools -t${nthreads} transform ${fullunion} histogram ${fullunion}.hist; cut -f2 ${fullunion}.hist | paste -sd+ | bc > ${cardloc}
  
  ## -ci1 -cs2
    
    '''
    def union_command(self,debug=False):
        pass

    
    def union_sketch(self, debug=False):
        ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmd = self.union_command()
        if not self.sketch_check():
            print(f"The sketch file {self._sfp.full} either doesn't exist or is empty")
            #self.cmd=subprocess.run(cmdlist)
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
        else:
            self.cmd = cmd
            #" ".join(cmdlist)
        if debug:
            print(self.cmd)
    ##TODO: make leaf sketch and union sketch private
        
    def create_sketch(self, debug=False):
        ''' If sketch file exists, assign path to self.sketch and return path. 
            If not create sketch, assign, and then return path.'''
        if self._sfp.ngen == 1:
            self.leaf_sketch(debug=debug)
        elif self._sfp.ngen > 1:
            self.union_sketch(debug=debug)
        else:
            raise RuntimeError("For some reason you are trying to sketch an empty list of files. Don't do that.")
        
        self.sketch = self._sfp.full
        return self.sketch
    

    def individual_card(self, speciesinfo:SpeciesSpecifics, debug:bool=False):
        '''Run cardinality for an individual sketch or database. Add it to a dictionary {path:value}'''
        pass
        
    def check_cardinality(self, speciesinfo: SpeciesSpecifics, debug=False):
        if self._sfp.full not in speciesinfo.cardkey.keys() or speciesinfo.cardkey[self._sfp.full] == 0:
            self.individual_card(speciesinfo=speciesinfo, debug=debug)   
        return float(speciesinfo.cardkey[self.sketch])

class DashSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[], delay=False, debug=True):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches, delay=delay, debug=debug)
    
    def individual_card(self, speciesinfo: SpeciesSpecifics, debug: bool = False):
        '''Run cardinality for an individual sketch. Add it to a dictionary {path:value}'''
        cmdlist = [DASHINGLOC,"card --presketched -p10"] +  [self.sketch]
        cmd = " ".join(cmdlist)
        if debug:
            print(cmd)
        #print(cmd)
        card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
        for card in csv.DictReader(card_lines, delimiter='\t'):
            speciesinfo.cardkey[card['#Path']] = card['Size (est.)']
        #return super().individual_card(speciesinfo, debug)

    def sketch_check(self):
        if (os.path.exists(self._sfp.full)) and os.stat(self._sfp.full).st_size == 0:
            return True
        else:
            return False
    
    def leaf_command(self):
        cmdlist = [DASHINGLOC, "sketch", 
        canon_command(self.experiment['canonicalize'], "dashing"),
        "-k" + str(self.kval), 
        "-S",str(self.experiment['registers']),
         "-p10","--prefix", str(self._sfp.dir),
          self._sfp.ffiles[0]]
                #    os.path.join(speciesinfo.inputdir, sfp.files[0])]
        cmd = " ".join(cmdlist)
        return cmd

    # def leaf_sketch(self, debug=False, delay=False, nthreads=10):
    #         ''' If leaf sketch file exists, record the command that would have been used. If not run the command and store it.'''
    #         cmdlist = [DASHINGLOC, "sketch", canon_command(self.experiment['canonicalize'], "dashing"),"-k" + str(self.kval),
    #         #TODO: get registers another way maybe
    #                 "-S",str(self.experiment['registers']),
    #                 "-p10","--prefix", str(self._sfp.dir),
    #                     self._sfp.ffiles[0]]
    #                 #    os.path.join(speciesinfo.inputdir, sfp.files[0])]
    #         cmd = " ".join(cmdlist)
    #         if debug:
    #                 print(cmd)
    #         if (not os.path.exists(self._sfp.full)) or os.stat(self._sfp.full).st_size == 0:
    #             print("The sketch file {0} either doesn't exist or is empty".format(self._sfp.full))
    #             subprocess.call(cmd, shell=True)
    #             self.cmd=cmd
    #             ##TODO Check if call raises error / returns 0
    #             ##TODO pass back the command so that it can be done in parallel?
    #         else:
    #             self.cmd = cmd
    def union_command(self, debug=False):
        cmdlist = [DASHINGLOC, "union", "-p 10 ","-z -o", str(self._sfp.full)] + self._presketches
        cmd = " ".join(cmdlist)
        return(cmd)
        
    def union_sketch(self, sfp, debug=False):
        ''' If union sketch file exists, record the command that would have been used. If not run the command and store it.'''
        cmdlist = [DASHINGLOC, "union", "-p 10 ","-z -o", str(sfp.full)] + self._presketches
        cmd = " ".join(cmdlist)
        #print("SIZE of {0} is {1}".format(self._sfp.full, os.stat(sfp.full)))
        if (not os.path.exists(sfp.full)) or os.stat(sfp.full).st_size == 0:
            # print("The sketch file {0} either doesn't exist or is empty".format(sfp.full))
            #self.cmd=subprocess.run(cmdlist)
            subprocess.call(cmd, shell=True)
            self.cmd=cmd
        else:
            self.cmd = cmd
            #" ".join(cmdlist)
        if debug:
            print(self.cmd)

class KMCSketchObj(SketchObj):
    def __init__(self, kval, sfp, speciesinfo, experiment, presketches=[], delay=False, debug=True):
        super().__init__(kval=kval, sfp=sfp, speciesinfo=speciesinfo, experiment=experiment, presketches=presketches, delay=delay, debug=debug)

    def individual_card(self, speciesinfo: SpeciesSpecifics, debug: bool = False):
        cmdlist = ["kmc_tools","info"] +  [self.sketch]
        cmd = " ".join(cmdlist)
        if debug:
            print(cmd)
        card_lines=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True).stdout.readlines()
        for line in card_lines:
            print(line)
            key, value = line.strip().split(':')
            print (f"{key},{value}")
            if key.strip() == 'total k-mers':
                speciesinfo.cardkey[self.sketch] = value.strip()
        #return super().individual_card(speciesinfo, debug)
    
    def sketch_check(self):
        if (os.path.exists(self._sfp.full +".kmc_pre")) and (os.path.exists(self._sfp.full +".kmc_suf")) and os.stat(self._sfp.full +".kmc_suf").st_size == 0:
            return True
        else:
            return False

    def leaf_command(self, nthreads=10, debug=False):
        tmpdir="/tmp/dandD"
        os.makedirs(tmpdir,exist_ok=True)
        cmdlist = ['kmc -t'+ str(nthreads),
        '-ci1 -cs2',f'-k{str(self.kval)}',
        canon_command(canon=self.experiment['canonicalize'], tool="kmc"),
        '-fm', self._sfp.ffiles[0], self._sfp.full, tmpdir]
        cmd = " ".join(cmdlist)
        return cmd

    # def leaf_sketch(self, delay=False, nthreads=10, debug=False):
    #     ''' If leaf kmc database exists, record the command that would have been used. If not run the command and store it.'''
    #     ##TODO: is this the right working directory to create/use?
    #     tmpdir="/tmp/dandD"
    #     os.makedirs(tmpdir,exist_ok=True)
    #     cmdlist = ['kmc -t'+ str(nthreads),'-ci1 -cs2',f'-k{str(self.kval)}',canon_command(canon=self.experiment['canonicalize'], tool="kmc"),'-fm', self._sfp.ffiles[0], self._sfp.full, tmpdir]
    #     cmd = " ".join(cmdlist)
    #     if debug:
    #             print(cmd)
    #     if (not os.path.exists(self._sfp.full +".kmc_pre")) or (not os.path.exists(self._sfp.full +".kmc_suf")) or os.stat(self._sfp.full +".kmc_suf").st_size == 0:
    #         print("The sketch file {0} either doesn't exist or is empty".format(self._sfp.full +".kmc_suf"))
    #         subprocess.call(cmd, shell=True)
    #         self.cmd=cmd

    def union_command(self, debug=False, nthreads=10):
        complex_input = "INPUT: \n"
        inputn=1
        for sketch in self._presketches:
            complex_input = complex_input + f"input{inputn} = {sketch} -ci1   \n"
            inputn+=1
        complex_input = complex_input + f"OUTPUT:\n{self._sfp.full} = " + " + ".join([f"input{i+1}" for i in range(inputn-1)])
        cmdlist = [f'echo -e "{complex_input}"',"|","kmc_tools","-t"+ str(nthreads), "complex", "/dev/stdin"]
        cmd = " ".join(cmdlist)
        return cmd
