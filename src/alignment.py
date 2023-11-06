import re
import numpy as np

class Alignment:
    def __init__(self, fname=None, alignment=None, root="naive") -> None:
        pass
        if fname is not None:
            alignment = self.read_fasta(fname)
            self.alignment= {key: list(value.strip()) for key,value in alignment.items()}
        else:
            self.alignment = alignment

        self.root= root
        self.reduced_pos =[]
        self.same_pos = []
        self.simplified=False
        self.full_root_seq = self.alignment[self.root]
        self.length = len(self.full_root_seq)
        


    @staticmethod
    def read_fasta(fname):
        '''read a fasta file into a dictionary 
        https://coding4medicine.com/backup/Python/reading-fasta-files.html '''
        seq = {}
        f=open(fname,'r')
        lines=f.readlines()

        hre=re.compile('>(\S+)')
        lre=re.compile('^(\S+)$')


        for line in lines:
                outh = hre.search(line)
                if outh:
                        id=outh.group(1)
                else:
                        outl=lre.search(line)
                        if(id in seq.keys()):
                                seq[id] += outl.group(1)
                        else:
                                seq[id]  =outl.group(1)
        return seq

    
    def write_fasta(self,fname):
            with open(fname, 'w') as file:
                    for key in self.alignment:
                            file.write(f">{key}\n")
                            file.write(f"{self.alignment[key]}\n")
    

    
    def simplify(self):

        self.simplified=True
        align_mat = np.array([self.alignment[key] for key in self.alignment])

        for j in range(align_mat.shape[1]):
            arr = align_mat[:,j]
            if np.unique(arr).shape[0] > 1:
                self.reduced_pos.append(j)
        
        # self.same_pos = [i for i in range(align_mat.shape[1]) if i not in self.reduced_pos ]
            
        align_mat = align_mat[:,self.reduced_pos]
        alignment = {}
        for i,key in enumerate(self.alignment.keys()):
              alignment[key] = align_mat[i,:].tolist()
        
        self.alignment = alignment
        return self.alignment
              


            
    def complete(self, mydict):

        alignment = {key: [] for key in mydict}
        pointer = 0
        for j in range(self.length):
            if j in self.reduced_pos:
                for seq in mydict:
                    alignment[seq].append(mydict[seq][pointer])
                pointer +=1
            else: 
                for seq in mydict:
                    alignment[seq].append(self.full_root_seq[j])
        

        labels = {key : "".join(value) for key, value in alignment.items()}
        return labels 
    
                     
                
                    

                
                
                
                
                
        

    


    
