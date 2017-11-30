Alignment_Algo <- function(seq_a, seq_b) {
 
 seq_a = "TCACACTA";
 seq_b = "AGCACACA";
 use = matrix(0,9,9);
 colnames(use) <- c("",substr(seq_a, 1, 1),substr(seq_a, 2, 2),substr(seq_a, 3, 3),substr(seq_a, 4, 4),substr(seq_a, 5, 5),substr(seq_a, 6, 6),substr(seq_a, 7, 7),substr(seq_a, 8, 8))
 rownames(use) <- c("",substr(seq_b, 1, 1),substr(seq_b, 2, 2),substr(seq_b, 3, 3),substr(seq_b, 4, 4),substr(seq_b, 5, 5),substr(seq_b, 6, 6),substr(seq_b, 7, 7),substr(seq_b, 8, 8))
 directions = use ;
#i =2;j=2;
 for (i in 2:9){ 
   for (j in 2:9){
     
     if (substr(seq_a, j-1, j-1) == substr(seq_b, i-1, i-1)){
       mat = use[i-1,j-1]+3
   }
     
     else {
       mat = use[i-1,j-1]-1
   }
     
     col_steady = use[i-1,j]-2
     row_steady = use[i,j-1]-2
     value = max(c(mat,col_steady,row_steady,0));
     use[i,j] = value;
     
     if (value == mat){
       directions[i,j]=3
     }
     
     else if (value == col_steady){   #back arrow goes up
       directions[i,j]=2          
     }
     
     else if (value == row_steady){   #back arrow goes left
       directions[i,j]=1
     }
   } 
     
 }  
 elem = use[9,9];
 score = elem;
 new_a="";
 new_b="";
 i = 9;
 j = 9;
 
 while(elem!=0) {
   
  arrow = directions[i,j];
  if (arrow == 3){
    new_a = paste(substr(seq_a, j-1, j-1),new_a,sep="");
    new_b = paste(substr(seq_b, i-1, i-1),new_b,sep="");
    #print(new_a);
    #print(new_b);
    i = i-1;
    j = j-1;
  }
  
  else if (arrow == 2){
    new_a = paste("_",new_a,sep="");
    new_b = paste(substr(seq_b, i-1, i-1),new_b,sep="");
    i = i-1;
  }
  
  else if (arrow == 1){
    new_a = paste(substr(seq_a, j-1, j-1),new_a,sep="");
    new_b = paste("_",new_b,sep="");
    j = j-1;
  }
  elem = use[i,j];
  
 }
 print(score);
 print(new_a);
 print(new_b);
}