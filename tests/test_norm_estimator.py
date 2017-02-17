import random
import numpy

random.seed(8675309)#set seed for reproducible results!
d= 10000;  #%big dimension
m=100;   #%small dimension
p = 120;  #%p is dimension of subspace.  m < p < d
#d = 36868752 #the number of junctions in tcga; it's big!
#m = 3000 #we have default dim 3000 for our morna indexes
#p = 10000 #an estimate for the size of "interesting" subspace
 
#%create m by d feature hashing matrix Phi
#Phi = zeros(m,d);
##phi = [ [0 for _ in range(d)] for _ in range(m)] #want m rows and d columns
phi = numpy.matrix([0]*(d*m))
phi.shape = (m,d)
#for j=1:d

for j in range(d):
#    I = randperm(m);
#    Phi(I(1),j) = sign(rand-.5);
	spot=random.randint(0,m-1) #randint is inclusive; avoid index out of bounds!
	##phi[spot][j] = random.choice((-1,1))
	phi[spot,j] = random.choice((-1,1))
#end


#%create d-dimensional test vector (normalized so that norm(x,2) = 1)
#x = zeros(d,1);
#x(1:p,1)=1;
#x(p+1:d,1) = .1;

x = numpy.matrix([1.0 for _ in range(d)] )
x = x.transpose()
for idx in range(p,d):
	x[idx,0]=.1
#x = x/norm(x,2);

x=x/numpy.linalg.norm(x,2)
 
#% y is feature hashing of x.  norm(y,2) should be rougly equal to
#% norm(x,2)
#y = Phi*x;
y = numpy.matmul(phi,x)
print("norm of y is " + str(numpy.linalg.norm(y,2)))
print("norm of x is " + str(numpy.linalg.norm(x,2)))

#norm(y,2);
 
 
#%e1 is a biased estimator for norm(x(1:p),2);
#%I believe that E[e1^2] = ((m-1)/m) norm(x(1:p,2),2)^2 + p/m
#z = Phi(:,1:d)'*y; #'* aka transpose (') and then matmult (*) each of the columns of phi with y
#is equivalent to dot product for column vector for each of these
z = numpy.matrix([0.0 for _ in range(d)] )
for column_num in range(d):
	z[0,column_num] = (phi[:,column_num].transpose() * y )
#e1 = norm(z(1:p,1),2)
e1 = numpy.linalg.norm(z[0,range(0,p)],2)
 
#%Thus using that norm(x,2)^2 = 1, 
#%I suggest exploring e2 = sqrt( max{0,e1^2 - p/m}) 
#%as a nearly unbiased
#%proxy for norm(x(1:p),2)
 
#if e1^2 > p/m
#    e2 = sqrt((m/(m-1))*(e1^2 - p/m));
#else
#    e2 = 0;
#end
e2 = sqrt( max(0,e1**2 - float(p)/m))
 
#%Compare 
#norm(x(1:p),2)  
#e2
print("True Norm: " + str(numpy.linalg.norm(x[range(p)])))
print("e2 Estimate: " + str(e2))