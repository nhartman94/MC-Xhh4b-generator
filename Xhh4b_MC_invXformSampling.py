import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ROOT import TLorentzVector
import sys
import argparse
import random
import os

'''
This function uses the pdf sin(x)/2, computed the cdf, which is a uniform random number in
the range [0,1], and returns the x such that cdf(x) = y.  

Note: For this generator, x is the polar angle, theta, in the rest frame of the decaying p'cle


'''
def pdf_sin():

	y = random.uniform(0.0, 1.0)
	x = np.arccos(1.0 - 2.0 * y)

	#  Here are some print statements to make sure that I get the same result whether I 
	#  use the analytic expression for x or the bisection method to find the root.
	#print "Analytic x = %1.2f" % x
	#
	## Since the cdf for sin(x) is cdf(x) = 1/2 * (1 - cos(x)), and as 
	## cdf(0) = 0 and cdf(pi) = 1, this means for a z in [0,1], cdf(x) - z
	## must have a zero intersected by the bisection method.
	#
	#z = y
	#f = lambda x: 0.5 * (1.0 - np.cos(x)) - z
	#print "Bisection x = %1.2f\n" % opt.bisect(f,0.0,np.pi)

	return x


'''
Use the pdf for theta for a spin 2 p'cle to generate a theta for a p'cle decaying in its rest frame.
'''
def pdf_sin4():

	z = random.uniform(0.0,1.0)
	
	# This is the function that we need to find the zero of (using the bisection method to) for sin^4(x) pdf
	f = lambda x: (x - 2.0/3.0*np.sin(2.0*x) + np.sin(4.0*x)/12.0) / np.pi - z

	return opt.bisect(f,0.0,np.pi)

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
Sorts two p'cles, pcle1 and pcle2 by their pT. 
The first p'cle in the output represents the p'cle with the larger pT.
'''
def pTOrder(pcle1,pcle2):

	if pcle1.Pt() > pcle2.Pt():
		return pcle1, pcle2
	else:
		return pcle2, pcle1


'''
This should be a way to simplify all the information that I'm getting ;) 
Make a subplot with 2 rows and 2 cols to plot pT, eta, phi, and m for a given input p'cle

Inputs: 	
	arr: an np array with N rows and 4 cols, where the rows are for the decay simulations,
	and the columns are the variables pT, eta, phi, and m

	Nbins: # of bins used to plot the histograms

	filepath: The filepath dictating the kinematics of the X decay simulated
	
	pcle: a string representing the particle for which plots are being made

This function plots the 4 histograms in 4 subplots, and saves the resultant figure using the 
filepath and var information.

'''
def plotKinematics(arr,Nbins,filepath,pcle,hist_range):

	fig, axes = plt.subplots(2,2)
	variables = ['$p_T$ (GeV)','$\eta$','$\phi$','m (GeV)']

	for i,ax in enumerate(fig.axes):

		# Label the plt
		if i == 0:
			ax.set_title('Distributions for ' + pcle)

		ax.hist(arr[:,i],Nbins,hist_range[i])
		ax.set_xlabel(variables[i])
		ax.set_ylabel('frequency')

		if hist_range[i] != None:
			ax.set_xlim(hist_range[i])
			if arr[:,i].min() < hist_range[i][0]:
				print("Warning: min({2} {3}) = {0} < hist_min = {1}".format(arr[:,i].min(), hist_range[i][0], pcle, variables[i]))
			if arr[:,i].max() > hist_range[i][1]:
				print("Warning: max({2} {3}) = {0} > hist_max = {1}".format(arr[:,i].max(), hist_range[i][1], pcle, variables[i]))

	plt.subplots_adjust(wspace=0.5, hspace=0.5)

	# Show the plt 
	#plt.show()

	# Save the plot
	fig.savefig( filepath + pcle + '.png' )
	plt.close()

	return True


'''
For a given Higgs, this function plots the dR b/w the leading and subleading higgs, the pT ratio b/w the leading
and subleading b, and a 2d hist of the dR b/w the bs vs. the pT of the h

Inputs:
	dR: The delta R between the leading and subleading b quarks
	pT_ratio: the ratio of the pT of the leading b to the subleading b
	pT_h: a numpy array containing the pT of the Higgs
	Nbins: the bins for which we're plotting the histograms
	filepath: the filepath to save the files in
	h_num: the number of the Higgs we're considering, where 1 is leading, 2 is subleading
'''
def plotAngularInfo(dR,pT_ratio,pT_h,Nbins,filepath,h_num):

	#plt.figure(1+3*(h_num-1))
	plt.hist(dR,Nbins)
	plt.title("$\Delta$R between the leading and subleading b-quarks")
	plt.xlabel("$\Delta$R(b{0}1,b{0}2)".format(h_num))
	plt.ylabel("Frequency")
	plt.savefig(filepath + "dR_h{}.png".format(h_num))
	plt.close()

	#plt.figure(2+3*(h_num-1))
	plt.hist(pT_ratio,Nbins,(0,100)) 
	plt.title("$p_T$ ratio between the leading and subleading b-quarks")
	plt.xlabel("$p_T$(b{0}1)/$p_T$(b{0}2)".format(h_num))
	plt.ylabel("Frequency")
	plt.savefig(filepath + "pT_ratio__h{}.png".format(h_num))
	plt.close()

	#plt.figure(3+3*(h_num-1))
	plt.hist2d(pT_h,dR,bins=Nbins,norm=LogNorm()) # Default bins is Nbins_x = Nbins_y = 10
	plt.colorbar()
	plt.xlabel("$p_T$(h{}) (GeV)".format(h_num))
	plt.ylabel("$\Delta$R(b{0}1,b{0}2)".format(h_num))
	plt.savefig(filepath + "dR_vs_pT_h{}.png".format(h_num))
	plt.close()

	return True


'''
This is the function that is called when the program is executed

It decays a scalar graviton in its rest frame into 

Inputs:
	mX:	mass of the graviton
	N:	The number of decay simulations run for this mass point
	pX:	momentum of the graviton 
'''
def run(mX, N, pX):

	print "mX =  %1.0f" % mX
	print "N iterations = %d" % N

	# Initialize 3 numpy arrays to hold the deltaR's that we're computing
	# and another 2 to hold the ratio b/w the leading and subleading bs in each h
	dR_hh = np.zeros(N)
	dR_h1 = np.zeros(N)
	dR_h2 = np.zeros(N)
	pT_ratio_h1 = np.zeros(N)
	pT_ratio_h2 = np.zeros(N)


	# The _arr suffix means this var is an Nx4 matrix, where the rows represent the decays,
	# and the columns stand for the kinematic variables: pT, eta, phi, and m
	X_arr   = np.zeros((N,4))  
        h1_arr  = np.zeros((N,4))
        h2_arr  = np.zeros((N,4))
        b11_arr = np.zeros((N,4))
        b12_arr = np.zeros((N,4))
        b21_arr = np.zeros((N,4))
        b22_arr = np.zeros((N,4))

	# Give the graviton some momentum in the eta = phi = 0 dir, i.e, along +x 
	X = TLorentzVector( pX, 0.0, 0.0, np.sqrt( np.square(pX) + np.square(mX) ) )
	vect_X = X.BoostVector()

	# Masses of the Higgs and b-quarks
	mh = 125.0
	mb = 5.0

	# Inline function for returning the kinematic vars for the 4-vector of a p'cle
	kinematics = lambda pcle: (pcle.Pt(), pcle.Eta(), pcle.Phi(), pcle.M())

	# Loop over the generated decays 
	for i in xrange(N):

		X_arr[i,0],X_arr[i,1],X_arr[i,2],X_arr[i,3] = kinematics(X)

		# Decay the graviton in its rest frame to two daughter particles (Higg's with masses 125 GeV)
		h1, h2 = decay(mX,mh,pdf_sin4)

		# Boost the h1, h2 back into the lab frame
		h1.Boost(vect_X)
		h2.Boost(vect_X)
		h1,h2 = pTOrder(h1,h2)		
		
		# Now that we have the 4-vector for the Higgses, use this to get the directions
		vect_h1 = h1.BoostVector()
		vect_h2 = h2.BoostVector()

		dR_hh[i] = h1.DeltaR(h2)
		
		h1_arr[i,0],h1_arr[i,1],h1_arr[i,2],h1_arr[i,3] = kinematics(h1)
		h2_arr[i,0],h2_arr[i,1],h2_arr[i,2],h2_arr[i,3] = kinematics(h2)

		# Then decay h1, h2 into b quarks (in the rest frame of the individual Higgses)
		b11, b12 = decay(mh,mb,pdf_sin)
		b21, b22 = decay(mh,mb,pdf_sin)

		# Boost the b-quark decay products into the lab frame
		b11.Boost(vect_h1)
		b12.Boost(vect_h1)
		b21.Boost(vect_h2)
		b22.Boost(vect_h2)

		# Sort the b-quarks by pT
		b11, b12 = pTOrder(b11,b12)
		b21, b22 = pTOrder(b21,b22)

		dR_h1[i] = b11.DeltaR(b12)
		dR_h2[i] = b21.DeltaR(b22)

		pT_ratio_h1[i] = b11.Pt() / b12.Pt()
		pT_ratio_h2[i] = b21.Pt() / b22.Pt()

		b11_arr[i,0],b11_arr[i,1],b11_arr[i,2],b11_arr[i,3] = kinematics(b11)
		b12_arr[i,0],b12_arr[i,1],b12_arr[i,2],b12_arr[i,3] = kinematics(b12)
		b21_arr[i,0],b21_arr[i,1],b21_arr[i,2],b21_arr[i,3] = kinematics(b21)
		b22_arr[i,0],b22_arr[i,1],b22_arr[i,2],b22_arr[i,3] = kinematics(b22)

	
	if(pX == 0):
		filepath = 'Xrest_mX_%dGeV_Nit_%d/' % (mX,N)
	else:
		filepath = 'pX_%dGeV_mX_%dGeV_Nit_%d/' % (pX,mX,N)

	if not os.path.isdir( os.path.join(os.getcwd(), filepath) ):
		os.mkdir(filepath)
	else:
		print("Warning: will overwrite the files in " + filepath) 


	# Make plots of the resulting kinematics
	Nbins = 100
	X_range = [(0,500),(-5.0,5.0),(-np.pi,np.pi),(0,1.2*mX)] 
	h_range = [None,   (-5.0,5.0),(-np.pi,np.pi),(0,1.2*mh)] 
	b_range = [None,   (-5.0,5.0),(-np.pi,np.pi),(0,1.2*mb)] 
	plotKinematics(X_arr,  Nbins,filepath,'X',  X_range)
	plotKinematics(h1_arr, Nbins,filepath,'h1', h_range)
	plotKinematics(h2_arr, Nbins,filepath,'h2', h_range)
	plotKinematics(b11_arr,Nbins,filepath,'b11',b_range)
	plotKinematics(b12_arr,Nbins,filepath,'b12',b_range)
	plotKinematics(b21_arr,Nbins,filepath,'b21',b_range)
	plotKinematics(b22_arr,Nbins,filepath,'b22',b_range)

	plotAngularInfo(dR_h1,pT_ratio_h1,h1_arr[:,0],Nbins,filepath,1)
	plotAngularInfo(dR_h2,pT_ratio_h2,h2_arr[:,0],Nbins,filepath,2)

	return True


# Read in the argument
parser = argparse.ArgumentParser()
parser.add_argument('--mX', help="float, the mass of the graviton", type=float)
parser.add_argument('--N',  help="int, # of decays simulated for each mass point", type=int)
parser.add_argument('--pX', help="float, momentum of the graviton, for now in the eta = phi = 0 dir", type=float)

args = parser.parse_args()

mX = args.mX if args.mX is not None else 500.0
N  = args.N  if args.N  is not None else 100
pX = args.pX if args.pX is not None else 0

run(mX,N,pX)

#for mX in [500, 1000, 1500, 2000]:
#	for pX in [0, 100, 200, 300]:
#		run(mX,N,pX)

