"""
 Names: Kevin Davis, Ackeem Greene & Semeion Stafford 
 Course: CMPS 4559
 Date: Nov/25/2019
 Description: This program Searches through the Mycoplasma genitalium
 DNA sequence for Open Reading Frames. It uses ATG as start codons and
 TAA, TAG, and TGA as stop codons. When a gene is found it prints its
 start index and its length. The ORF translation will only be printed
 for the first gene found in every offset.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable
table = CodonTable.ambiguous_dna_by_id[1]

output= open("output.txt","w+")
print(" Names: Kevin Davis, Ackeem Greene & Semeion Stafford \n",
      "Course: CMPS 4559\n",
      "Date: Nov/25/2019\n",
      "Description: This program Searches through the Mycoplasma genitalium\n",
      "DNA sequence for Open Reading Frames. It uses ATG as start codons and\n",
      "TAA, TAG, and TGA as stop codons. When a gene is found it prints its\n",
      "start index and its length. The ORF translation will only be printed\n",
      "for the first gene found in every offset.",
      file = output)

for seq_record in SeqIO.parse("mycoplasma.gb","genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

# use this loop to switch to reverse complement
for x in range (2):
    #take seq_record and store it in mySec
    mySec = (seq_record.seq )#,IUPAC.unambiguous_dna)
    
    #on the second iteration, we will search through the reverse complement
    if x > 0:
        mySec = mySec.reverse_complement()
        print("\n\n\n\nReverse Complement\n",
              "#######################################################\n",
              "#######################################################", 
              file = output)
        
    gene_info =[]
    #print out everything in mySec, character by character
    for index, letter in enumerate(mySec):
        gene_info.append(letter)
    found_start = False
    found_end = False
    start_index = -1
    start_end_distance = 0
    gene = ""
    counter = 0
    offSet = 0
    codon = [0,0,0]

    # this loop cycles through the 3 offsets
    while offSet < 3:
        #used to only print out the translation first ORF in sequence
        First_ORF = 0
        print ("\n\n\noffSet: ", offSet,
               "\n-------------------------------------------------------\n",
               file = output)
        
        #this loop cycles through to find the start and stop codons
        # -3 + offset to not go out of range
        while counter * 3 < len(gene_info) - (3 + offSet):
            codon[0] = gene_info[counter * 3 + + offSet]
            codon[1] = gene_info[counter * 3 + 1 + offSet]
            codon[2] = gene_info[counter * 3 + 2 + offSet]
            #store codons in a string
            codon_str = str(codon[0] + codon[1] + codon[2])
            #if no start codon has been found yet
            if found_start == False and codon_str == "ATG" :
                found_start = True
                #add one because counter start from 0 but document start from 1
                start_index = counter * 3 + 1 + offSet     
            #if found, found start will be true, increment distance
            if found_start == True:
                #distance starts at 1, not 0
                start_end_distance = start_end_distance + 3
                #concatenate codons together to make the gene
                gene = gene + codon_str
            #if an end is found but the distance is less than 1000, start over
            if found_start == True and start_end_distance < 1000 and (codon_str == "TAG" or codon_str == "TAA" or codon_str == "TGA"):
                #reset values
                found_end = False
                found_start = False
                start_end_distance = 0
                start_index = -1
                gene = ""
            #if an end is found and the dstance is 1000 or more  
            if found_start == True and start_end_distance >= 1000 and (codon_str == "TAG" or codon_str == "TAA" or codon_str == "TGA"):
                 #output gene information
                 #only print out the translation for the first ORF found in every offset
                 # I assume this is what you wanted from reading the guidelines
                 if First_ORF == 0:
                     #convert the string of gene to a sequence
                     SeqGene = Seq(gene,generic_dna)
                     new_gene = SeqGene.translate(table = "Bacterial", to_stop = True)
                     #tring the translation
                     print (new_gene, file = output)
                 print ("Start Index is :", start_index, file = output)
                 print ("Distance From Start :", start_end_distance, "\n", file = output)
                 #reset values to find another gene
                 found_end = False
                 found_start = False
                 start_end_distance = 0
                 start_index = -1
                 gene = ""
                 #if a gene is found increment First ORF
                 First_ORF += 1
            counter += 1
            
        #when checking another offset reset all values to default
        found_end = False
        found_start = False
        start_end_distance = 0
        start_index = -1
        counter = 0
        gene = ""
        #increment offset
        offSet +=1
    
    
output.close()