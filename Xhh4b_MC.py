import numpy as np
import scipy
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
This is the function where I'm going to put the decay function for the spin2 graviton
'''
def spin2():

	return 1.0


'''
The decay function that describes how a particle decays
Inputs:
	m_parent:  	mass of the parent
	m_daughter:	mass of the daughters
	dist:   The distribution function for the decay of p'cle X
	^ I'll include the distribution dependence later


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
	
	phi, theta = dist()

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

	# Initialize the graviton to start at rest (for now)
	X = TLorentzVector(0.0, 0.0, 0.0, mX)
	
	# Masses of the Higgs and b-quarks
	mh = 125.0
	mb = 5.0

	# Put a for loop here
	for i in range(N):

		# Decay the graviton in its rest frame to two daughter particles (Higg's with masses 125 GeV)
		h1, h2 = decay(mX,mh,uniform)

		# Now that we have the 4-vector for the Higgses, use this to get the directions
		vect_h1 = h1.BoostVector()
		vect_h2 = h2.BoostVector()

		dR_hh[i] = h1.DeltaR(h2)
		
		# We decay in the rest frame of mX, and then we can boost it (or h1,h2) into the lab frame
		# Note, at this stage we're not including the spin of the parent, or the fact that it's not at rest... 

		# Then decay h1, h2 into b quarks
		b11, b12 = decay(mh,mb,uniform)
		b21, b22 = decay(mh,mb,uniform)

		# Boost the b-quark decay products into the lab frame
		b11.Boost(vect_h1)
		b12.Boost(vect_h1)
		b21.Boost(vect_h2)
		b22.Boost(vect_h2)

		
		dR_h1[i] = b11.DeltaR(b12)
		dR_h2[i] = b21.DeltaR(b22)


	# Make plots of the resulting kinematics
	# Goal: look at how the b-quarks become collimated as you vary the mass of X

	plt.hist(dR_hh,N/10)
	plt.title("$\Delta$R between the two Higgses")
	plt.xlabel("$\Delta$R")
	plt.ylabel("Frequency")
	plt.show()
	

	plt.hist(dR_h1,N/10)
	plt.title("$\Delta$R between the b-quarks from the 1st Higgs")
	plt.xlabel("$\Delta$R")
	plt.ylabel("Frequency")
	plt.show()

	plt.hist(dR_h2,N/10)
	plt.title("$\Delta$R between the b-quarks from the 2nd Higgs")
	plt.xlabel("$\Delta$R")
	plt.ylabel("Frequency")
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
	run(mX_default,N_default)
