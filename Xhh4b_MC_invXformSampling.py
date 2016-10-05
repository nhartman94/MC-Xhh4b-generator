import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from ROOT import TLorentzVector
import sys
import argparse
import random


'''
This function outputs phi and theta taken from a distribution uniform in cos(theta) and phi

'''
def uniform():

	# Assume uniform in (x,y,z) => uniform in d(cos(theta)), phi
	phi    = random.uniform(0.0, 2.0*np.pi) 
	cos_th = random.uniform(-1.0,1.0) 

	# From cos_th, you can use the arccos to get back to theta
	theta = np.arccos(cos_th)

	return phi, theta


'''
This function uses the pdf sin(x)/2, computed the cdf, which is a uniform random number in
the range [0,1], and returns the x such that cdf(x) = y.  

Note: For this generator, x is the polar angle theta in the rest frame of the decaying p'cle
'''
def pdf_sin():

	y = random.uniform(0.0, 1.0)
	x = np.arccos(1.0 - 2.0 * y)

	return x


'''
Use the pdf for theta for a spin 2 p'cle to generate a theta for a p'cle decaying in its rest frame.


'''
def pdf_sin4():

	global z  # Need to modify global copy of the global variable z
	z = random.uniform(0.0,1.0)
	return opt.bisect(f,0.0,np.pi)

'''
This is the function that we need to find the zero of (using the bisection method to) for stuff
'''
def f(x):
	return (x - 2.0/3.0*np.sin(2.0*x) + np.sin(4.0*x)/12.0) / np.pi - z


'''
The decay function that describes how a particle decays
Inputs:
	m_parent:  	mass of the parent
	m_daughter:	mass of the daughters
	dist:   The pdf for theta for the decay of p'cle X

Outputs: the TLorentz 4-vectors for the two daughters
	v1: 	the daughter p'cle with higher pT
	v2: 	the daughter p'cle with lower  pT


'''
def decay(m_parent, m_daughter, dist):

	# Somehow need to use the dist to describe how the daughter particles decay...
	
	# Since this function assumes that the parent is at rest, the energy for each of the daughter is half of the parent's mass
	E = m_parent / 2.0

	# Initialize the daughter TLorenz vectors to 0s.
	v1 = TLorentzVector(1.0, 0.0, 0.0, E)


	# All decays considered in this generator are symmetric about the z-axis => they'll be uniform in phi
	phi   = random.uniform(0.0, 2.0*np.pi) 
	theta = dist()

	v1.SetRho( np.sqrt( E**2 - m_daughter**2))
	v1.SetTheta(theta)
	v1.SetPhi(phi)
	
	# Define v2 using Cartesian coordinates, i.e, px,py,pz 
	v2 = TLorentzVector(-v1.Px(),-v1.Py(),-v1.Pz(),E)

	return v1, v2



'''
This is the function that is called when the program is executed

It decays a scalar graviton in its rest frame into 

Inputs:
	mx:	mass of the graviton
	N:	The number of decay simulations run for this mass point

'''
def run(mX, N):

	print "mX =  %1.0f" % mX
	print "N iterations = %d" % N

	# Initialize 3 numpy arrays to hold the deltaR's that we're computing
	dR_hh = np.zeros(N)
	dR_h1 = np.zeros(N)
	dR_h2 = np.zeros(N)


	# Also initialize numpy arrays for costh and phi for the variables.
	costh_h1  = np.zeros(N) 
	phi_h1    = np.zeros(N)  
	costh_b11 = np.zeros(N) 
	phi_b11   = np.zeros(N)  
	costh_b12 = np.zeros(N) 
	phi_b12   = np.zeros(N)  
	costh_h2  = np.zeros(N) 
	phi_h2    = np.zeros(N)  
	costh_b21 = np.zeros(N) 
	phi_b21   = np.zeros(N)  
	costh_b22 = np.zeros(N) 
	phi_b22   = np.zeros(N)  

	# Initialize the graviton to start at rest (for now)
	X = TLorentzVector(0.0, 0.0, 0.0, mX)
	
	# Masses of the Higgs and b-quarks
	mh = 125.0
	mb = 5.0

	# Put a for loop here
	for i in range(N):

		# Decay the graviton in its rest frame to two daughter particles (Higg's with masses 125 GeV)
		#h1, h2 = decay(mX,mh,uniform)
		h1, h2 = decay(mX,mh,pdf_sin4)

		# Now that we have the 4-vector for the Higgses, use this to get the directions
		vect_h1 = h1.BoostVector()
		vect_h2 = h2.BoostVector()

		dR_hh[i] = h1.DeltaR(h2)
		costh_h1[i] = h1.CosTheta()
		phi_h1[i]   = h1.Phi()
		costh_h2[i] = h2.CosTheta()
		phi_h2[i]   = h2.Phi()

	
		# We decay in the rest frame of mX, and then we can boost it (or h1,h2) into the lab frame
		# Note, at this stage we're not including the spin of the parent, or the fact that it's not at rest... 

		# Then decay h1, h2 into b quarks
		b11, b12 = decay(mh,mb,pdf_sin)
		b21, b22 = decay(mh,mb,pdf_sin)

		# Boost the b-quark decay products into the lab frame
		b11.Boost(vect_h1)
		b12.Boost(vect_h1)
		b21.Boost(vect_h2)
		b22.Boost(vect_h2)


		dR_h1[i] = b11.DeltaR(b12)
		dR_h2[i] = b21.DeltaR(b22)


		costh_b11[i] = b11.CosTheta()
		phi_b11[i]   = b11.Phi()
		costh_b12[i] = b12.CosTheta()
		phi_b12[i]   = b12.Phi()
		costh_b21[i] = b21.CosTheta()
		phi_b21[i]   = b21.Phi()
		costh_b22[i] = b22.CosTheta()
		phi_b22[i]   = b22.Phi()


	# Make plots of the resulting kinematics
	# Goal: look at how the b-quarks become collimated as you vary the mass of X

	#plt.hist(dR_hh,N/10)
	#plt.title("$\Delta$R between the two Higgses")
	#plt.xlabel("$\Delta$R")
	#plt.ylabel("Frequency")
	##plt.show()
	#

	#plt.hist(dR_h1,N/10)
	#plt.title("$\Delta$R between the b-quarks from the 1st Higgs")
	#plt.xlabel("$\Delta$R")
	#plt.ylabel("Frequency")
	##plt.show()

	#plt.hist(dR_h2,N/10)
	#plt.title("$\Delta$R between the b-quarks from the 2nd Higgs")
	#plt.xlabel("$\Delta$R")
	#plt.ylabel("Frequency")
	##plt.show()

	# Plot the distribution in phi and cos_th for h1, h2, and b11, b12, b21, b22
	# Plot b11, b22 on the same histogram => watch them become collimated??
	Nbins = N/50

	plt.hist(costh_h1, Nbins,alpha=1,facecolor='w',edgecolor='k',label='h1')
	plt.hist(costh_h2, Nbins,alpha=1,facecolor='w',edgecolor='b',label='h2')
	plt.legend(loc='upper right')
	plt.xlabel('cos($\Theta$)')
	plt.ylabel('Frequency')
	plt.title('The polar angle distribution for the Higgses')
	plt.show()

	plt.hist(costh_b11,Nbins,alpha=1,facecolor='w',edgecolor='r',label='b11')
	plt.hist(costh_b12,Nbins,alpha=1,facecolor='w',edgecolor='g',label='b12')
	plt.legend(loc='upper right')
	plt.xlabel('cos($\Theta$)')
	plt.ylabel('Frequency')
	plt.title('The polar angle distribution from the b\'s from the \"leading\" h')
	plt.show()


	plt.hist(costh_b21,Nbins,alpha=1,facecolor='w',edgecolor='c',label='b21')
	plt.hist(costh_b22,Nbins,alpha=1,facecolor='w',edgecolor='m',label='b22')
	plt.legend(loc='upper right')
	plt.xlabel('cos($\Theta$)')
	plt.ylabel('Frequency')
	plt.title('The polar angle distribution from the b\'s from the \"subleading\" h')
	plt.show()

	plt.hist(phi_h1, Nbins,alpha=1,facecolor='w',edgecolor='k',label='h1')
	plt.hist(phi_h2, Nbins,alpha=1,facecolor='w',edgecolor='b',label='h2')
	plt.hist(phi_b11,Nbins,alpha=1,facecolor='w',edgecolor='r',label='b11')
	plt.hist(phi_b12,Nbins,alpha=1,facecolor='w',edgecolor='g',label='b12')
	plt.hist(phi_b21,Nbins,alpha=1,facecolor='w',edgecolor='c',label='b21')
	plt.hist(phi_b22,Nbins,alpha=1,facecolor='w',edgecolor='m',label='b22')
	plt.legend(loc='upper right')
	plt.xlabel('$\phi$')
	plt.ylabel('Frequency')
	plt.show()

	return True


# Read in the argument
parser = argparse.ArgumentParser()
parser.add_argument('--mX', help="float, the mass of the graviton", type=float)
parser.add_argument('--N', help="int, # of decays simulated for each mass point", type=int)

args = parser.parse_args()

mX_default = 500.0
N_default = 100

# Pass arguments to main
if args.mX is not None and args.N is not None:
        run(args.mX, args.N)
elif args.mX is not None:
	run(args.mX,N_default)
elif args.N is not None:
	run(mX_default,args.N)
else:
	# plot the variables as a function of mass points (?)
	run(mX_default,N_default)
