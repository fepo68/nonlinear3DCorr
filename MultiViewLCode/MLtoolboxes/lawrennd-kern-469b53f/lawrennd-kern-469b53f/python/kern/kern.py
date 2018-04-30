import sys
import ndlutil
from ndlutil.utilities import sigmoid
import numpy as np
import pdb

class kern(ndlutil.parameterised):
	"""The base class for kernels. Should not be instantiated"""
	def __init__(self,X):
		ndlutil.parameterised.__init__(self)
		self.shape = (X.shape[0],X.shape[0])
		self.dimensions = X.shape[1]
		self.masked=False
		self.scaled=False
		self.name='unnamed'
	def set_param(self):
		raise NotImplementedError
	def get_param(self):
		raise NotImplementedError
	def get_param_names(self):
		raise NotImplementedError
	def set_X(self):
		raise NotImplementedError
	def cross_args(self,X2):
		"The arguments to the function when computing the covariance between self.X and X2"
		raise NotImplementedError
	def gradients(self):
		raise NotImplementedError
	def expand_X(self,X):
		"""apply masks and scaling to X and call set_X"""
		self.Xoriginal = X.copy()
		self.shape=(X.shape[0],X.shape[0])
		if self.scaled:
			X = self.scalefunc(X)
		if self.masked:
			self.set_X(X[self.mask,:])
		else:
			self.set_X(X)
		
	def set_mask(self,mask):
		"""Apply a mask to this kernel. 

		Arguments
		---------
		mask : a valid index that can be applied to self.X: a boolean array of sefl.X.shape[0] or an integer array
		"""
		if mask.size == self.X.shape[0]:
			mask = np.nonzero(mask.flatten())[0]
		if mask.size==self.X.shape[0]:
			self.masked=False
			return
		self.masked=True
		self.mask = mask
		self.mask_grid = np.meshgrid(mask,mask)
		self.expand_X(self.Xoriginal)
	
	def set_scaled(self,limits=None):
		"""force the kernel to scale the data before computing.
		convenintly deals with scaling args for cross computing etx.
		
		Arguments
		---------
		limits = a tuple of (xmin, xmax)
		"""
		if limits is None:
			self.Xmin,self.Xmax = self.Xoriginal.min(), self.Xoriginal.max()
		else:
			self.Xmin,self.Xmax = limits
		self.scaled=True
		self.scalefunc = lambda x : (x-self.Xmin)/(self.Xmax-self.Xmin)
		self.expand_X(self.X)

	def compute(self,target=None):
		"""
		Arguments
		---------
		target : a np array to add the computation of this kernel to, for in-place computation. If None, a new array is created

		Notes
		-----
		This takes care of the masking of the kernel, as well as in-place computation.
		"""
		if target is None:
			target = np.zeros(self.shape)
		if self.masked:
			target[self.mask_grid] += self.function(*self.args)
		else:
			target += self.function(*self.args)
		return target

	def cross_compute(self,X2,target=None):
		"""Compute the covariance between self.X and the passed values X2. Mask aware."""
		if target is None:
			target = np.zeros((self.shape[0],X2.shape[0]))
		if self.scaled:
			args = self.cross_args(self.scalefunc(X2))
		else:
			args = self.cross_args(X2)
		if self.masked:
			target[self.mask,:] += self.function(*args)
		else:
			target += self.function(*args)
		return target

	def compute_new(self,Xnew):
		"""Compute the covariance of an unrelated set of points, Xnew. """
		X = self.X
		shape = self.shape
		masked = self.masked
		if masked:
			self.masked=False
		self.expand_X(Xnew)
		ret = self.compute()
		if masked:
			self.masked=True
		self.set_X(X)
		self.shape = shape
		return ret
		
	def diag_compute(self):
		"""Compute just the diagonal terms of the covariance"""
		raise NotImplementedError, "TODO"


	def extract_gradients_new(self,targets=None):
		x = self.get_param()
		gradients = self.gradients(*self.args)

		#work out which gradient go where, with overlaps due to tying and nans due to fixing
		target_index = np.arange(x.size)
		for tie in self.tied_indices:
			target_index[tie[1:]] = tie[0]
		some_fixed=False
		if len(self.constrained_fixed_indices):
			some_fixed=True
			target_index[np.hstack(self.constrained_fixed_indices)] = target_index.size

		target_index = np.argmax(target_index[:,None]==np.unique(target_index),1) #NB: the fixed indices are now just the maximum value

		#get the gradient correction factors which some from the constraints TODO: arbitrary constraints, not just exp
		factors = np.ones(x.size)
		[np.put(factors,i,x[i]) for i in self.constrained_positive_indices]
		[np.put(factors,i,x[i]) for i in self.constrained_negative_indices]
		[[np.put(factors,i,(x[i]-l)*(h-x[i])/(h-l)) for i in inds] for inds,h,l in zip(self.constrained_bounded_indices,self.constrained_bounded_uppers, self.constrained_bounded_lowers)]

		if targets is None:
			targets = [np.zeros(self.shape) for i in range(np.unique(target_index).size-some_fixed)]
		if self.masked:
			for f,g,i in zip(factors,gradients,target_index):
				if (not some_fixed)|(i<target_index.max()):
					targets[i][self.mask_grid] += f*g
		else:
			[np.add(targets[ti],g*f,targets[ti]) for ti,g,f in zip(target_index,gradients,factors) if (not some_fixed)|(ti<target_index.max())]
		return targets
		


	def extract_gradients(self,targets=None):
		"""
		Extract the gradients. Much like model.extract_gradients,
		but we're dealing with lists of matrices here, which is a little more tricky.

		Also deals with the masking of the kernel
		
		Also does in-place computation!

		Returns
		-------
		A list of np arrays.
		"""
		x = self.get_param()
		if targets is None:
			targets = [np.zeros(self.shape) for i in range(x.size)]
		if self.masked:
			#[np.add(t[self.mask_grid],g,t[self.mask_grid]) for t,g in zip(targets,self.gradients(*self.args))]
			for t,g in zip(targets,self.gradients(*self.args)):
				t[self.mask_grid] += g
		else:
			[np.add(t,g,t) for t,g in zip(targets,self.gradients(*self.args))]


		x = self.get_param()
		#in-place multiplication correct gradients. in-place summation for tied parameters. 
		[np.multiply(targets[i],x[i],targets[i]) for i in self.constrained_positive_indices]
		[np.multiply(targets[i],x[i],targets[i]) for i in self.constrained_negative_indices]
		[[np.multiply(targets[i],(x[i]-l)*(h-x[i])/(h-l),targets[i]) for i in inds] for inds,h,l in zip(self.constrained_bounded_indices,self.constrained_bounded_uppers, self.constrained_bounded_lowers)]
		[[np.add(targets[ii], targets[i[0]], targets[i[0]]) for ii in i[1:]] for i in self.tied_indices]
		to_remove = self.constrained_fixed_indices + [t[1:] for t in self.tied_indices]
		if len(to_remove):
			to_remove = np.hstack(to_remove)
		else:
			to_remove = []
		return [gg for i,gg in enumerate(targets) if not i in to_remove]

	def extract_gradients_X(self,target=None):
		"""
		The gradient of the kernel wrt the input points self.X
		wraps the implementation of gradients_X(), utilizing the mask

		Notes
		-----
		In-place computation available by passing the target array. 

		Returns
		-------
		A NxNxD array. the i,j,kth element is the gradient of the [i,j]the element of the kernel wrt the k'th dimension of the i'th input.
		"""
		if target is None:
			target = np.zeros(self.shape+(self.dimensions,))
		if self.masked:
			#this is pretty wierd: why do the indeces swap like that?
			target[self.mask_grid[1],self.mask_grid[0],:] += self.gradients_X()
		else:
			target += self.gradients_X()
		return target

	def __add__(self,other):
		if isinstance(other,compound):
			other.__add__(self)
		else:
			return compound(self.X,[self,other])
class compound(kern):
	def __init__(self,X,kerns=[]):
		kern.__init__(self,X)
		if isinstance(kerns,kern):
			kerns = [kern]

		self.kerns = []
		self.shape = kerns[0].shape
		self.param_allocation = []
		[self.__add__(k) for k in kerns]
		self.args = tuple()#needed so that we can utilise extract_gradients

	def set_param(self,x):
		[k.expand_param(x[i]) for k,i in zip(self.kerns, self.param_allocation)]
	def get_param(self):
		return np.hstack([k.extract_param() for k in self.kerns]+[np.zeros(0)])#this does a safe hstack
	def get_param_names(self):
		return sum([['cmpnd%i_%s_%s'%(i,k.name,n) for n in k.extract_param_names()] for i,k in enumerate(self.kerns)],[])

	def set_X(self,X):
		[k.expand_X(X) for k in self.kerns]

	def compute(self,target=None):
		if target is None:
			target = np.zeros(self.shape)
		[k.compute(target) for k in self.kerns]
		return target

	def cross_compute(self,X2,target=None):
		if target is None:
			target = np.zeros((self.shape[0],X2.shape[0]))
		[k.cross_compute(X2,target) for k in self.kerns]
		return target

	def compute_new(self,Xnew):
		return np.sum([k.compute_new(Xnew) for k in self.kerns],0)

	def gradients(self):
		return sum([k.extract_gradients() for k in self.kerns],[]) #simply concatenates all the gradient matrices!

	def gradients_X(self):
		ret = np.zeros(self.shape+(self.kerns[0].X.shape[1],))
		[k.extract_gradients_X(ret) for k in self.kerns] # in place computation
		return ret

	def __add__(self,other):
		if isinstance(other,compound):
			assert self.shape==other.shape
			self.kerns.extend(other.kerns)
			Nparam = self.get_param().size
			self.param_allocation.extend([i+Nparam for i in other.param_allocation])
			print 'Warning! parameter constraints may be lost! (TODO)'# TODO
		elif isinstance(other,kern):
			Nparam = self.get_param().size
			self.kerns.append(other)
			assert self.shape==other.shape
			self.param_allocation.append(np.arange(Nparam,Nparam+other.extract_param().size))
		else:
			raise AttributeError, "%s is not a proper kernel"%str(other)
		return self

class cross_kern(kern):
	def __init__(self,kern1, kern2):
		assert (kern1.masked and kern2.masked), "mask the crossed kernels before instantiating the cross kern object"
		kern.__init__(self, kern1.Xoriginal)
		self.kern1, self.kern2 = kern1, kern2
		self.masked=True
		self.name = self.kern1.name +'X' + self.kern2.name
		#self.maskgrid = (self.kern # TODO!!
		#TODO: deal with scaling?
	def set_X(self,X1, X2):
		"""
		Any precomputation on X can be done here by children of the class
		"""
		raise NotImplementedError
	def cross_args(self,X2):
		"The arguments to the function when computing the covariance between self.X and X2"
		raise NotImplementedError # This needs some thought
	def expand_X(self,X):
		"""
		Take the matrix X and the meshes from the crossed kernels and precompute the with set_X
		"""
		# Take the masks from the linked kernels.
		self.set_X(X[self.kern1.mask,],X[self.kern2.mask]) 
		
	def set_mask(self):
		"""
		Take the mask from the crossed kernels
		"""
		self.maskgrid = (self.kern1.maskgrid[0],self.kern2.maskgrid[1])

	def compute(self):
		if target is None:
			target = np.zeros(self.shape)
		f = self.function(*self.args)
		target[self.mask_grid] += f
		target[self.mask_grid[::-1]] += f # no transopse: reversing the mesh does this.
		return target


class structured(compound):
		def __init__(self,X,kerns,conn, xkerns='all'):
			"""
			A kernel with structured covariance 
			This is more general than hierarhcical, below, but it _can_ do hierarhcical structures as well as multi-ouptut structures

			Arguments
			---------
			X - (array, N x D), the continuos domain of the kernel 
			kerns - (list of kern, len Dout) the kernels, each on X
			conn - (N X Dout binary array) the connectivity of the kernels. 
			xkerns - (Dout x Dout binary array), or 'None', or 'all':  which cross kernels to construct
			
			"""

			#make sure things are the right size
			self.Dout = len(kerns)
			assert connections.shape[1] == self.Dout
			assert xkerns in ['all', None] or xkerns.shape == (self.Dout,self.Dout)

			self.connectivity = conn
			self.X = X

			#initialize and mask the kernels
			compound.__init__(self,X,kerns)
			for i,k in enumerate(self.kerns):
				assert k.shape==kerns[0].shape
				k.set_mask(np.nonzero(connections[:,i])[0])

			#sort out the cross kernels
			if xkerns == 'all':
				xkerns = np.ones((self.Dout,self.Dout))
			elif xkerns is None:
				xkerns = np.zeros((self.Dout, self.Dout))
			assert np.all(xkerns==xkerns.T) # make sure xkern connectivity is symmetrical
			for i,j in itertoos.combinations(range(self.Dout,2)):
				if xkerns[i,j]:
					self.kerns.append(find_cross_kerns(self.kerns[i],self.kerns[j]))

		def cross_compute(self,X2,conn2):
			raise NotImplementedError, "TODO"

		def compute_new(self,Xnew, conn):
			raise NotImplementedError, "TODO"

		def __add__(self,foo):
			raise AttributeError, "Adding kernels to a structured kern is not currently supported"
			

class hierarchical(compound):
	def __init__(self,X,connections,kerns):
		"""
		Arguments
		---------
		X : the numpy array on which the kernel is to be computed
		connections : a numpy boolean array with as many rows as X.
		kerns : a list of kernels  with the same domain X.

		Notes
		-----
		the n,k th element of connections is True if the nth input is affected by the kth kernel

		Example
		-------
		con = np.hstack((np.ones((90,1)),np.kron(np.eye(3),np.ones((30,1)))))
		X = np.linspace(0,10,90)[:,None]
		kerns = [kern.rbf(X) for i in range(4)]
		k = kern.hierarchical(con,kerns)

		"""
		assert len(kerns)==connections.shape[1]
		compound.__init__(self,X,kerns)
		self.connections = connections
		for i,k in enumerate(self.kerns):
			assert k.shape==kerns[0].shape
			mask = np.nonzero(connections[:,i])[0]
			if len(mask)<connections.shape[0]:
				k.set_mask(mask)

	def cross_compute(self,X2,connections):
		if len(self.kerns)==1:
			connections = connections.reshape(connections.size,1)
		assert connections.shape == (X2.shape[0],len(self.kerns)), "Connections do not match input data"
		ret = np.zeros((self.shape[0],X2.shape[0]))
		for k,c in zip(self.kerns,connections.T):
			i = np.nonzero(c)[0]
			if len(i):
				ret[:,i] += k.cross_compute(X2[i,:]) #can't do in place here due to fancy indexing (which always copies) TODO: in-place would save some time
		return ret

	def compute_new(self,Xnew,connections):
		if len(self.kerns)==1:
			connections = connections.reshape(connections.size,1)
		assert connections.shape == (Xnew.shape[0],len(self.kerns)), "Connections do not match input data"
		ret = np.zeros((Xnew.shape[0],Xnew.shape[0]))
		for k,c in zip(self.kerns,connections.T):
			i = np.nonzero(c)[0]
			if len(i):
				ret[np.meshgrid(i,i)] += k.compute_new(Xnew[i,:]) 
		return ret
			
	
	def __add__(self,kerns,connections=None):
		"""
		Add another kernel to this structure. 

		Fully connected by default
		"""
		
		if isinstance(kerns,kern):
			kerns = [kerns]
		for i,k in enumerate(kerns):
			compound.__add__(self,k)
			if not connections is None:
				k.set_mask(np.nonzero(connections[:,i])[0])
			

	




