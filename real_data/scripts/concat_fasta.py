import argparse
import re



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

def write_fasta(fname, mydict):
        with open(fname, 'w') as file:
                for key in mydict:
                        file.write(f">{key}\n")
                        file.write(f"{mydict[key]}\n")

if __name__ =="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-l", "--light")
    parser.add_argument("-y", "--heavy")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    light_seqs = read_fasta(args.light)
    heavy_seqs = read_fasta(args.heavy)

    assert set(light_seqs.keys()) == set(heavy_seqs.keys())
    concat_seqs = {}
    for key, value in heavy_seqs.items():
        concat_seqs[key] = value + light_seqs[key]
    
    sorted_dict = {key: value for key, value in sorted(concat_seqs.items())}
    
    write_fasta(args.outfile, sorted_dict)





