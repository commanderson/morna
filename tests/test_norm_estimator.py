import argparse
import numpy
import random
import sys
from math import sqrt


#First set up parser and parse in input file's name.
parser = argparse.ArgumentParser()
parser.add_argument('-s', metavar='<int>', type=int,
            required=False,
            default=1,
            help='The number of blocks the hashing space is broken into'
        )
        
parser.add_argument('-l', metavar='<int>', type=int,
            required=False,
            default=10,
            help='The number of re-drawing loops of the projection matrix'
        )
        
args = parser.parse_args()


random.seed(8675309)#set seed for reproducible results!
d= 10000;  #%big dimension
m=100;   #%small dimension
p = 120;  #%p is dimension of subspace.  m < p < d
#Note that the following settings are WAY TOO MUCH to handle:
#d = 36868752 #the number of junctions in tcga; it's big!
#m = 3000 #we have default dim 3000 for our morna indexes
#p = 10000 #an estimate for the size of "interesting" subspace
 


#%create d-dimensional test vector (normalized so that norm(x,2) = 1)
#x = zeros(d,1);
#x(1:p,1)=1;
#x(p+1:d,1) = .1;

#Old vector generation code exactly following matlab follows:
#x = numpy.matrix([1.0 for _ in range(d)] )
#x = x.transpose()
#for idx in range(p,d):
#	x[idx,0]=.1
#x = x/norm(x,2);

#New vector generation of pseudo-random vector
x = numpy.matrix([random.uniform(1.0,1000.0) for _ in range(d)] )
x = x.transpose()

#normalize x for norm of 1
x=x/numpy.linalg.norm(x,2)
 
real_values = []
estimations = []
e1s = []
xnorm = []
ynorm = []
blocksize = m/args.s

for loop in range(args.l):
	if not (loop % 100):
		sys.stdout.write("Loop %d\r" % loop)
	#%create m by d feature hashing matrix Phi
	#Phi = zeros(m,d);
	##phi = [ [0 for _ in range(d)] for _ in range(m)]
	#want m rows and d columns
	phi = numpy.matrix([0.0]*(d*m))
	phi.shape = (m,d)
	#for j=1:d

	for j in range(d):
		for block in range(args.s):
			#print("Block %d is from %d to %d" %
			#		 (block,block*blocksize,((block+1)*blocksize)-1))
			spot=random.randint(block*blocksize,((block+1)*blocksize)-1) 
			#randint is inclusive; avoid index out of bounds!
			##phi[spot][j] = random.choice((-1,1))
			phi[spot,j] = random.choice((-1.0/sqrt(args.s),1.0/sqrt(args.s)))
			print("just set [phi%d, %d] to %f" % (spot,j,phi[spot,j]))
			#phi[spot,j] = random.choice((-1,1))
	#end

	#% y is feature hashing of x.  norm(y,2) should be rougly equal to
	#% norm(x,2)
	#y = Phi*x;
	y = phi * x
	#y = numpy.matmul(phi,x)
	#y=y/sqrt(args.s)
	#print("norm of y is " + str(numpy.linalg.norm(y,2)))
	#ynorm.append(numpy.linalg.norm(y,2)/sqrt(args.s))
	ynorm.append(numpy.linalg.norm(y,2))
	#print("norm of x is " + str(numpy.linalg.norm(x,2)))
	xnorm.append(numpy.linalg.norm(x,2))
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
 	#e1 = e1/sqrt(args.s)
 	e1s.append(e1)
	#%Thus using that norm(x,2)^2 = 1, 
	#%I suggest exploring e2 = sqrt( max{0,e1^2 - p/m}) 
	#%as a nearly unbiased
	#%proxy for norm(x(1:p),2)
 
	#if e1^2 > p/m
	#    e2 = sqrt((m/(m-1))*(e1^2 - p/m));
	#else
	#    e2 = 0;
	#end

	#e2 = sqrt( max(0,(float(m)/(m-1))*(e1**2 - float(p)/m)))
	#print(e1**2 - float(p)/m)
	e2 = sqrt( max(0,(e1**2 - float(p)/m)))
	#e2 = e2/sqrt(args.s)
	#%Compare 
	#norm(x(1:p),2)  
	#e2


	#print("True Norm: " + str(numpy.linalg.norm(x[range(p)])))
	real_values.append(numpy.linalg.norm(x[range(p)]))
	#print("e2 Estimate: " + str(e2))
	estimations.append(e2)
	#print("Ratio TN/e2: " + str(numpy.linalg.norm(x[range(p)])/e2))

sys.stdout.write("\n")
differences = [real_values[i]-estimations[i] for i,_ in enumerate(estimations)]	
e2_rv_ratio = [estimations[i]/real_values[i] for i,_ in enumerate(estimations)]	
#print estimations
print("mean real value:" + str(numpy.mean(real_values)))
#print("real value variance:" + str(numpy.var(real_values)))
print("mean estimation:" + str(numpy.mean(estimations)))
print("estimation variance:" + str(numpy.var(estimations)))
print("mean difference:" + str(numpy.mean(differences)))
print("difference variance:" + str(numpy.var(differences)))
print("mean e2/rv ratio:" + str(numpy.mean(e2_rv_ratio)))

print("xnorm:" + str(numpy.mean(xnorm)))
print("mean ynorm:" + str(numpy.mean(ynorm)))
print("ynorm variance:" + str(numpy.var(ynorm)))

#print("mean e1:" + str(numpy.mean(e1s)))
#print("e1 variance:" + str(numpy.var(e1s)))