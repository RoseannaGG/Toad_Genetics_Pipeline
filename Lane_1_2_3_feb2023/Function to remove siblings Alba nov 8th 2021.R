
namesib9<-read.table("F:/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/805_populationsruns_nov2021/805_populations_R70pctoverall_p19_minmaf0.01_maxhet0.6_writesinglesnp/601INDIV_602SNPS/Rmsibs/listindiv_kingrobust_0.25_r1_1.9.txt")

namesib9$V1

listsibs<-namesib9$V1


### as a for loop
for(i in 1:length(namesib9)){
     
      gl.toad_test <- gl.toad_test[indNames(gl.toad_test) != listsibs[i]]
       
      }
gl.toad_test@ind.names


## as a function

remove_sibs <- function(genlight_object,
                      names){
  
  initial_length  <- length(genlight_object@ind.names)
  
  print(paste0("Initial object made up of ", initial_length, " individuals"))
  
  for(i in 1:length(names)){
    
    genlight_object <- genlight_object[indNames(genlight_object) != names[i]]
    
  }  
  
 
  
  
  final_length <- length(genlight_object@ind.names)
  
  print(paste0("Final object made up of ", final_length, " individuals"))
  
  
  return(genlight_object)
  
  
}



### run it on data with function
#test_funk <- remove_sibs(genlight_object = gl.toad,
                        # names = listsibs)
