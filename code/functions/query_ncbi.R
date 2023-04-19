query_ncbi <- function(marker, organisms){
     # Given a region of interest ('marker') formulated as query term, searches
     # NCBI's Nucleotide database for all sequences corresponding to that query 
     # combined with all organisms listed in 'organisms'.  
     
     # Returns a DNAStringSet of results with full header preserved for 
     # any downstream manipulation/renaming. 
     
     # Requires the following packages: BioStrings, rentrez
     library(Biostrings)
     library(rentrez)
     
     # Long queries cause HTTP failure, so chunk into 50-species units
     organisms <- split(organisms, ceiling(seq_along(organisms)/50))
     
     organism_query <- 
          organisms %>%
          sapply(dQuote) %>%
          sapply(paste0, '[ORGN]') %>%
          sapply(paste, collapse = ' OR ') %>%
          sapply(function(x){paste0('(', x, ')')})

     # Add marker
     query <- paste(marker, organism_query, sep = ' AND ')
     
     ids <- lapply(query, 
                   entrez_search, 
                   db='nucleotide', retmax = 10000, use_history=TRUE)
     
     # Maybe-- flag how many times we hit retmax?
     # What's the default retmax?
     # Could also just prioritize trying this on the "old" version of human-foods
     
     seqs.return = NULL # Store results here
     seqs.count = 0
     
     for (i in seq_along(ids)){
          if (ids[[i]]$count==0){
               # Do nothing further if there aren't any hits
          } else {
               # Regular expression to pull FASTA header out from string
               # Don't match >, match immediately after up to first space
               ex <- '[^>]\\S*' 
               seqs <- 
                    entrez_fetch(db='nucleotide', 
                                 web_history = ids[[i]]$web_history, 
                                 rettype='fasta') %>%
                    # This returns concatenated sequence strings; split apart 
                    # so we can re-name inline
                    strsplit('\n{2,}') %>% # Usually two newline chars, but sometimes more
                    unlist()
               # Save this to ultimately combine with taxonomy data, as want to
               # be able to identify these sequences after the fact
               accs <- str_extract(seqs, ex) 
               
               # Keep full header for descriptive name
               headers <- str_extract(seqs, '^[^\n]*')
               
               seqs <- 
                    seqs %>%
                    # Now update seqs to sequence only, stripping header
                    sub('^[^\n]*\n', '', .) %>%
                    # And removing separating \n characters
                    gsub('\n', '', .)
               
               # Now add to DNAStringSet
               seqs <- DNAStringSet(seqs)
               names(seqs) <- headers
               
               # Update counters
               seqs.count <- seqs.count + ids[[i]]$count
               seqs.return <- append(seqs.return, seqs)
               
               cat('50 species processed... \n')
               cat('Equal lengths:', seqs.count == length(seqs.return), '\n')
          }
     } 
     
     cat(seqs.count, 'sequences processed for', marker, '\n')
     seqs.return
}