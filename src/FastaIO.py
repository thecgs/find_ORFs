#!/usr/bin/env python
# coding: utf-8

import re
import os
import io
import bz2
import sys
import gzip
import zipfile
import tarfile
import warnings
from CodonTable import getCodonTable

__author__ = "Guisen Chen <thecgs001@foxmil.com>"
__all__ = ['get_file_obj', 'get_files_obj', 'FastaIO', 'FastaSeqence', 
           'BaseSeqence', 'NucleicSeqence', 'ProteinSeqence','guess_fasta_type']

def get_file_obj(in_file):
    """Return a file object from an input file.
    """
    if not os.path.exists(in_file) and in_file != "-":
        raise Exception("can't open {}".format(in_file))
    elif in_file == '-':
        return sys.stdin
    if in_file.find(".tar") > 0:
        if in_file.endswith(".tar.gz"):
            tp = tarfile.open(in_file, "r:gz")
        elif in_file.endswith(".tar"):
            tp = tarfile.open(in_file, "r")
        elif in_file.endswith(".tar.bz2"):
            tp = tarfile.open(in_file, "r:bz2")
        return io.TextIOWrapper(tp)
    elif in_file.endswith(".gz"):
        return gzip.open(in_file, "rt")
    elif in_file.endswith(".zip"):
        zobj = zipfile.ZipFile(in_file)
        zp = zobj.open(zobj.namelist()[0], "r")
        return io.TextIOWrapper(zp)
    elif in_file.endswith(".bz") or in_file.endswith(".bz2"):
        return bz2.BZ2File(in_file, "rt")
    else:
        return open(in_file,'r')
    
    
def get_files_obj(in_files):
    """Accept a files list object, and return a iters _io.TextIOWrapper list object.
    """
    return [get_file_obj(infile) for infile in in_files]
    
    
def guess_fasta_type(seqence):
    """判断fasta file的文件类型，氨基酸序列返回True, 核酸序列返回False
    """
    res = False
    if ("M" in seqence) or ("V" in seqence) or ("L" in seqence)     or ("I" in seqence) or ("P" in seqence) or ("F" in seqence)     or ("Y" in seqence) or ("W" in seqence) or ("S" in seqence)     or ("Q" in seqence) or ("D" in seqence) or ("E" in seqence)     or ("K" in seqence) or ("R" in seqence) or ("H" in seqence) :
        res = True
    return res
    
    
class FastaIO:
    """FastaIO parser
    """
    def __init__(self, files):
        if not isinstance(files, list):
            self._files = get_files_obj([files])
        else:
            self._files = get_files_obj(files)
        self._strat = -1
        
    def __len__(self):
        return len(self._files)
    
    def __str__(self):
        return str(self._files)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        self._strat += 1
        if len(self._files) <= self._strat:
            raise StopIteration
        else:
            return self._files[self._strat]
    
    def __getitem__(self, item):
        return self._files[item]
    
    def __enter__(self):
         return self
    
    def __exit__(self):
        for file in self._files:
            file.close()

    def __del__(self):
        del self
        return None
        
    def names(self):
        return [file.name for file in self._files]
    
    def parse(self):
        for file in self._files:
            for lineno, line in enumerate(file):
                if  line.strip().startswith('>'):
                    if lineno != 0:
                        if guess_fasta_type(Seq):
                            messege = '\033[91m' + Name +'\033[93m is a protein sequence. Please check your sequence and enter a nucleic acid sequence.\033[0m'
                            warnings.warn(messege, category=Warning)
                            #print(Name +' is a protein sequence. Please check your sequence and enter a nucleic acid sequence.')
                            #yield FastaSeqence((ID, Name, Description, ProteinSeqence(Seq)))
                        else:
                            yield FastaSeqence((ID, Name, Description, NucleicSeqence(Seq)))
                    Name = line.strip()[1:]
                    ID = line.strip().split()[0][1:]
                    Description = Name[len(ID)+1:]
                    Seq = ""
                elif line.strip() == "":
                    pass
                else:
                    Seq += line.strip()
            yield FastaSeqence((ID, Name, Description, NucleicSeqence(Seq)))

class FastaSeqence:
    """返回一个fasta file格式的序列类
    """
    def __init__(self, args):
        ID, Name, Description, Seq = args
        self.ID = ID
        self.Name = Name
        self.Description = Description
        self.Seq = Seq
        
    def __len__(self):
        return len(self.Seq)
    
    def __getitem__(self,item):
        return self.Seq[item]
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return "[ID: {}\n Name: {}\n Description: {}\n Seq: {}]".format(self.ID, self.Name, self.Description, self.Seq[:20]+"......" + self.Seq[-20:])
    

class BaseSeqence:
    """返回一个序列类，是所有序列类的基类
    """
    def __init__(self, string):
        self._Seq = string
        self._type = 'General'
        
    def __len__(self):
        return len(self._Seq)
    
    def __getitem__(self, item):
        return self._Seq[item]
    
    def count(self, string):
        return self._Seq.count(string)
    
    def __str__(self):
        return str(self._Seq)
    
    def __repr__(self):
        return str((self._Seq[:20] + "......" + self._Seq[-20:], self._type))


class NucleicSeqence(BaseSeqence):
    """返回一个核苷酸序列类
    """
    def __init__(self, string):
        self._RNA2DNA = str.maketrans('Uu', 'Tt')
        self._Seq = str(string).translate(self._RNA2DNA)
        self._type = 'Nucleic'
    
    def GC(self):
        return self.count('G') + self.count('C')
    
    def AT(self):
        return self.count('A') + self.count('T')
    
    def ATGC(self):
        return self.AT() + self.GC()
    
    def GC_skew(self):
        return round((self.count('G') - self.count('C'))/self.GC(), 4)
        
    def upper(self):
        return NucleicSeqence(self._Seq.upper())
    
    def lower(self):
        return NucleicSeqence(self._Seq.lower())
        
    def translate(self, codontable=1):
        codontable=getCodonTable(codontable)[0]
        protein = ''
        n = int(len(self._Seq)//3)*3
        for site in range(0, n, 3):
            if "N" in self._Seq[site:site+3]:
                protein += "X"
            else:
                protein += codontable[self._Seq.upper()[site:site+3]]
        return ProteinSeqence(protein)
    
    def reverse(self):
        return NucleicSeqence(self._Seq[::-1])
    
    def compliment(self):
        intab = "ATGCatgc"
        outtab = "TACGtacg"
        trantab = str.maketrans(intab, outtab)
        return NucleicSeqence(self._Seq.translate(trantab))
    
    def reverse_compliment(self):
        return NucleicSeqence(str(self.compliment())[::-1])
    
class ProteinSeqence(BaseSeqence):
    """返回一个氨基酸序列类
    """
    def __init__(self, string):
        self._Seq = string
        self._type = 'Protein'
        
    def remove_StopCodon(self):
        """移去终止密码子
        """
        if self._Seq.endswith('.') or self._Seq.endswith('*'):
            return ProteinSeqence(self._Seq[:-1])
        else:
            return ProteinSeqence(self._Seq)

