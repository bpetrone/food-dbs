read_fastas <- function(path){
     # Given multiple FASTA files present in a directory, read them in and 
     # append their sequences to a shared DNAStringSet object.`
     files <- list.files(path, pattern = '.fna|.fa|.fasta',
                         full.names = TRUE)
     all_seqs <- DNAStringSet()
     
     for (file in files){
          seqs <- readDNAStringSet(file)
          all_seqs <- append(all_seqs, seqs)
     }
     
     all_seqs
}