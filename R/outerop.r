#' @export

#OUTEROP calculate an outer operation on two vectors
#   function y=outerop(A,B,OPERATOR)
# 
#   Calculates resultant matrix when the OPERATOR is applied 
#   to all combinations of the elements of vector A and the 
#   elements of vector B e.g. the outer product of A and B 
#   is outerop(A,B,'*'), the outer sum of A and B 
#   is outerop(A,B,'+') 
#
#   If OPERATOR is omitted '+' is assumed
#
#   This function is equivalent to the 
#   APL language's circle.dot operation. e.g. 
#   in APL Ao.*B is the outer product 
#   of A and B % and Ao.+B is the outer sum.
#   Ideally it would work on matrices but in practice
#   I've only ever needed to use it with vectors.
#
#   Copyright Murphy O'Brien 2005
#   all rights unreserved
#
outerop<-function(a,b,operator)
{
if (missing(operator))
    {operator="+"}                      # for only two arguments assume outerproduct   

if (operator=="*")                      # common operator 
    {y=a%*%b} else    
  {outera=as.matrix(a)%*%rep(1,length(b))          # these are the matrices that 
  outerb=as.matrix(rep(1,length(a)))%*%b # meshgrid(A,B) would generate 
  functionHandle=match.fun(operator)    # new R14 functionality
  y=functionHandle(outera,outerb)  }    # allows faster/neater method
return (y)
}

