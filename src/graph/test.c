#include <element.h>
#include <hash_value.h>

int main(int argc, char **argv){
  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  } 
  
  //allocate memory for the sliding windows         

  windows->window = malloc(sizeof(KmerSlidingWindow) * max_windows);       
  if (windows->window== NULL){
    fputs("Out of memory trying to allocate an array of KmerSlidingWindow",stderr);
    exit(1);
  }
  
  //allocate memory for every every sliding window
  int w;
  for(w=0;w<max_windows;w++){
    KmerSlidingWindow * current_window =&(windows->window[w]);
    
    current_window->kmer = malloc(sizeof(BinaryKmer) * max_kmers);
    if (current_window->kmer == NULL){
      fputs("binary_kmer: Out of memory trying to allocate an array of BinaryKmer",stderr);
      exit(1);
    }      
  }      
  //----------------------------------
  

  
}
