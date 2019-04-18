"Corner" <- function(y,x,Nrow=11,Ncol=7){
### compute the "corner-table" for output "y" and input "x".
### y: filtered dependent variable
### x: filtered input variable AND is supposed to be a white noise series
if(is.matrix(y))y=y[,1]
if(is.matrix(x))x=x[,1]
nT <- min(length(y),length(x))
lag <- Nrow+Ncol+1
y1 <- y[1:nT]
x1 <- x[1:nT]
Sy <- sqrt(var(y1))
Sx <- sqrt(var(x1))
Y <- y1[lag:nT]
X <- x1[lag:nT]
for (i in 1:(lag-1)){
X <- cbind(X,x1[(lag-i):(nT-i)])
}
vi <- cor(Y,X)
vi <- vi*Sy/Sx
vmax <- max(abs(vi))
vi <- vi/vmax
##cat("Corr: ",vi,"\n")
tbl <- matrix(0,Nrow, Ncol)
tbl[,1] <- vi[1:Nrow]
for (j in 2:Ncol){
for (i in 1:Nrow){
 cmx = diag(rep(vi[i],j))
 for (ii in 2:j){
  for (jj in 1:(ii-1)){
   idx = i+jj
   cmx[ii,jj] = vi[idx]
   }
  }
 for (jj in 2:j){
  for (ii in 1:(jj-1)){
   idx = i-jj+1
   if(idx > 0)cmx[ii,jj]=vi[idx]
   }
  }
#  cat("cmx: ",cmx,"\n")  
tbl[i,j]=det(cmx)
 }
 }
stbl <- tbl
crit=2/sqrt(nT)
tbl=cbind(c(1:Nrow)-1,tbl)
colnames(tbl) <- c("r->",paste(c(1:Ncol)))
cat("Corner Table: ","\n")
print(round(tbl,3))
for (i in 1:Nrow){
 for (j in 1:Ncol){
  if(abs(tbl[i,j+1]) <= crit){
      stbl[i,j] = "O"
      }
      else{
      stbl[i,j] = "X"
      }
   }
 }
cat("\n")
cat("Simplified Table: 2/sqrt(T): ","\n")
J=paste(c(1:Nrow)-1)
stbl=cbind(J,stbl)
colnames(stbl) <- c("r->",paste(c(1:Ncol)))
print(stbl)

Corner <- list(cornor = tbl)
}