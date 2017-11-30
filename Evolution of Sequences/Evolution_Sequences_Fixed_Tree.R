
#  Simulate evolution of sequences on a fixed tree  #

#installed.packages()["ape","Package"]
#install.packages("ape")
library(ape)

# define the tree

?read.tree

mytree<-read.tree(text="(orangutan:13,(gorilla:10.25,(human:5.5,chimp:5.5):4.75):2.75);")

#plot(mytree)

summary(mytree)
str(mytree)
mytree$edge
mytree$tip.label
mytree$edge.length
mytree<-reorder(mytree,order="cladewise")

# Construct the instantaneous rate matrix under TN93 model

## Definition of stationary distribution and kappa (transition/transversion ratio)
kappa1=44.229
kappa2=21.781
piT=0.22
piC=0.26
piA=0.33
piG=0.19


Q = matrix(c(-(kappa1*piC+piA+piG),kappa1*piC,piA,piG,
             kappa1*piT,-(kappa1*piT+piA+piG),piA,piG,
             piT,piC,-(piT+piC+kappa2*piG),kappa2*piG,
             piT,piC,kappa2*piA,-(piT+piC+kappa2*piA)),nrow=4,ncol=4,byrow=TRUE)



# Scale/Normalize it, so that the overall mutation rate is 1 per site (per unit of time - defined by tree)

normconst<- -(piT*Q[1,1] + piC*Q[2,2] + piA*Q[3,3] + piG*Q[4,4])
Q<-Q*(1/normconst)

#Simulate evolution of 20 sites, using your rate matrix

## Define the three different overall mutation rates
# in units of #mutations per site per million of years (MY) = unit of time in the tree
scalesvector<-c(0.000135,0.0135,1.35)

## Define a structure storing each of the 3 outputs (for each low, intermediate, high mutation rates)
results<-rep(list(rep(list(NULL),4)),3)

## Define the loop for three overall mutation rates
for(scales in scalesvector){

  Qnew<- Q*scales
  nucleotidefreq<-c(piT,piC,piA,piG)
  set.seed(1)

  ## Simulate evolution starting from root to the tips for 20 sites
  sitenumber<-20
  for (site in 1:sitenumber){
    # Create a vector to hold states of each (internal and external) node
    sitepattern<-rep(0,times=length(mytree$edge[,1])+1)

    ## Select what would be the first nucleotide based on the stationary distribution
    u<-runif(n=1,min = 0,max = 1)
    # Use  encoding: 1=T,2=C,3=A,4=G
    nucleotide<-0
    # Variable "sumfreq" counts the sum of the nucleotide stationary frequencies
    sumfreq<-c(piT,piC+piT,piC+piT+piA,1)
    # Find a nucleotide for which sumfreq_previous < u <= sumfreq_next
    if (u<=sumfreq[1]){
      nucleotide = 1
    }else if (u<=sumfreq[2]){
      nucleotide = 2
    }else if (u<=sumfreq[3]){
      nucleotide = 3
    }else{
      nucleotide = 4
    }
    sitepattern[mytree$edge[1,1]]<-nucleotide

    # Simulate the sequence evolution noting down the state at each internal/external node

    for(j in 1:length(mytree$edge[,1])){
      ## Simulate sequence changing from root down to the tips
      brlength<-mytree$edge.length[j]
      nucleotide<-sitepattern[mytree$edge[j,1]]

      ## Keep evolving until we exceed the branch length
      while(brlength>0){

        # Choose exponentially distributed waiting time until a change happens

        utime<- -( 1 / (-Qnew[nucleotide,nucleotide]))* log(runif(n=1,min = 0,max = 1))

        # Adjust branch length remaining

        brlength<- brlength - utime;

        # Check whether a change happens

        if(brlength>0){

          t=0;

          for (k in 1:4) {
            if (k!=nucleotide){
              t=t+1
              sumfreq[t] = Qnew[nucleotide,k] / -Qnew[nucleotide,nucleotide];
            }
          }
          sumfreq[2] = sumfreq[2]+sumfreq[1];
          sumfreq[3] = 1;
          if (nucleotide==1){
            possible = c(2,3,4)
          }
          else if (nucleotide==2){
            possible = c(1,3,4)
          }
          else if (nucleotide == 3) {
            possible = c(1,2,4)
          }
          else{
            possible = c(1,2,3)
          }

          u<-runif(n=1,min = 0,max = 1)

          if (u<=sumfreq[1]){
            nucleotide = possible[1]
          }
          else if (u<=sumfreq[2]){
            nucleotide = possible[2]
          }
          else{
            nucleotide = possible[3]
          }
        }
      }

      ## Record the final nucleotide at each internal/external node
      sitepattern[mytree$edge[j,2]]<-nucleotide
    }

    ## Record the full pattern for the given site in the results list, but only for external nodes
    for (extnode in 1:length(mytree$tip.label))
    results[[which(scalesvector==scales)]][[extnode]]<-c(results[[which(scalesvector==scales)]][[extnode]],sitepattern[extnode])
  }

  for (extnode in 1:length(mytree$tip.label)){
    result<-c(as.character(results[[which(scalesvector==scales)]][[extnode]]))
    cat(result,mytree$tip.label[extnode], "\n",file = paste("alignment_",scales,".txt",sep=""), sep = " ", fill = FALSE, labels = NULL,
    append = FALSE)
  }

  # translate the data back to A/C/G/T
  translate<-function(x){ if (x==1) return('T'); if (x==2) return('C'); if (x==3) return('A'); if (x==4) return('G')}
  filename = paste("alignment_",scales,".txt",sep="")
  cat("species alignment \n",file = filename, sep = "", fill = TRUE)
  for (extnode in 1:length(mytree$tip.label)){
    result<-results[[which(scalesvector==scales)]][[extnode]]
    result<-paste(unlist(lapply(result,translate)),collapse="")
    cat(mytree$tip.label[extnode], result, "\n", file = filename, sep = " ", fill = FALSE, labels = NULL,
    append = TRUE)
  }
}
