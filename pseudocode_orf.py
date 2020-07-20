def parse_command_line_options():

    ```
    This function will return a command line option
    ```
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta')
    parser.add_argument('--gff3')
    parser.add_argument('--output')
    parser.add_argument('--min_length',type=int)
    parser.add_argument('--start_codons',default = ['ATG'],nargs="+")
    parser.add_argument('--stop_codons',default = ['TAA','TAG','TGA'],nargs="+")
    args = parser.parse_args()
    
    return args
    
    
def extract_UTR5_CDS_postion(mygff3):
   
    ``` 
    Each sequence’s in .fasta file start position of UTR5, end position of UTR5 and end position of CDS will be extracted and written in to position_dict dictionary.
    ```
    
    db = gffutils.create_db(mygff3, 'myGFF.db', merge_strategy="create_unique", keep_order=True, force=True)
    db = gffutils.FeatureDB('myGFF.db')
    position_dict = {}
    for i in db.features_of_type("UTR5"):
        dict = {}
        for j in db.features_of_type("CDS"):
            if i.seqid == j.seqid:
               dict[i.seqid]=[i.start-1,j.end,i.end]
               position_dict.update(dict)
    
    return position_dict


def extract_fasta_sequence(seq,position_dict)
        
        ```
        This function extracts the UTR5 and CDS region sequence in the input .fasta file and writes them into a new fasta file
        ```
        
        for key in position_dict.keys(): 
            for record in SeqIO.parse(myfasta, "fasta"):     
               if key == record.id:                 
                  sequences = SeqRecord(Seq(str(record.seq[position_dict[key][0]:position_dict[key][1]])),id = record.id)
                  with open(newfasta, "a") as output_handle:
                        SeqIO.write(sequences, output_handle, "fasta")
                        
                        
def read_new_fasta(newfasta):

    ```
    This function reads the new fasta file and generates a fasta file iterator
    ```
    
    fasta_iterator = SeqIO.parse(newfasta, "fasta")
    
    return fasta_iterator
                        
                                                
def find_startcodon():

    ```
    
    This function will scan the string one by one just like ribosome. 
    
    It first determines if current codon is one of the three start codons -- ATG GTG TTG. If it is, its position will be recorded. Then to find the stop codon for the start 
    
    codon. Once the function find the first stop codon, the inner loop breaks. In this way, it avoids the orfs with in-frame stop codons.
    
    ```
    
    for record in SeqIO.parse(newfasta,"fasta"):
         seq = record.seq
         for frame in range(3):
             for i in range(frame, len(seq), 3):
                  current_codon1 = seq[i:i+3]                            
                  if current_codon1 in ["ATG","GTG","TTG"]:
                      start_codon_index = i 
                      for j in range(start_codon_index,len(seq),3):
                          current_codon2 = seq[j:j+3]
                          if current_codon2 in ["TAA","TAG","TGA"]:
                            stop_codon_index = j
                            length = stop_codon_index - start_codon_index + 3
                            start_uorf = start_codon_index + 1 
                            stop_uorf = stop_codon_index+3
                            yield (record.id, frame, start_uorf, stop_uorf,length)
                            break 
    

    
def check_length():

    ```
    This fuction will check if the length of orf meets the users' min_length parameter
    ```
    val = [list(ele) for ele in list(startstop_codon())]
    condition = lambda x: x[4] >= min_length * 3
    filter_list = list(filter(condition, val))
    
    return filter_list
    

def check_uorf_type():
    
    ```
    This funtion will check the which type the  uorf  belongs to and write the type element into the list
    ```
    
       val = check_length(10)
       for key in position_dict.keys():         
           for i in val:
              if key == i[0]:             
                 end_UTR5 = position_dict[key][2]                
                 end_CDS = position_dict[key][1]               
                 start_uorf = i[2]
                 stop_uorf = i[3]
                 if start_uorf < end_UTR5 and stop_uorf < end_UTR5:
                      i.append("uORF")
                 if start_uorf < end_UTR5 and stop_uorf == end_CDS:
                      i.append("CDS_NTE")
                 if start_uorf < end_UTR5 and end_UTR5 < stop_uorf < end_CDS:
                      i.append("overlap_uORF")

       return val
    
    
def remove_main_orfs():

    ```
    This function will remove orfs in the list who are not uorfs
    ```
    
    condition = lambda x: len(x) == 6
    filter_list = list(filter(condition, a))
    filter_list
    
 
def write_to_gff3():

  ```
  This function will write the list into a gff3 file
  ```
  
   for record in SeqIO.parse(newfasta,"fasta"):
    for i in filter_list:
        if record.id == i[0]:
           out_file = "test.gff3"                
           seq = record.seq #do not have real meaning , just for create a object
           name = i[0]
           rec = SeqRecord(seq, name)
           qualifiers = {"source": "riboviz", 
                         "score": ".", 
                         "start_codon": seq[i[2]-1:i[2]+2], 
                         "Name": name + "_" + i[5] + "_" + str(i[2]),
                         "frame":i[1]}                
           feature = SeqFeature(FeatureLocation(i[2]-1, i[3]), 
                                type= i[5], 
                                strand=1,
                                qualifiers=qualifiers)
           rec.features = [feature]
           with open(out_file, "a") as out_handle:                     
                 GFF.write([rec], out_handle)
   
    os.system("sed -i '/#/d' {}" .format("test.gff3"))   


def main():

    ```
    This function takes care of other function
    ```
