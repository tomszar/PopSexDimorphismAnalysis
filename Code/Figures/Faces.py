class MultipleFaceShapes(object):
	landmark_matrix = 0
	facets  = 0 
	normals = 0
	sh      = 0
	n_land  = 0
	n_obs   = 0

#Try to create an object for multiple faceshapes

class FaceShape(object):
	landmarks = 0
	facets    = 0
	normals   = 0
	sh        = 0
	n_land    = 0
	color     = 0
	
	def __init__(self, landmarks, facets, normals = 0, sh = (7160, 3), n_land = 0, color = (0.7,0.7,0.7)):
		self.landmarks = landmarks
		self.facets    = facets
		self.normals   = normals
		self.sh        = sh
		self.n_land    = n_land
		self.color     = color
		self.check_matrices()
		
	def check_matrices(self):
		import numpy as np
		
		#To numpy array
		self.landmarks = np.array(self.landmarks)
		self.facets    = np.array(self.facets)
		
		#Checking shape of landmark matrix
		if 1 in self.landmarks.shape:
			self.landmarks = np.reshape(self.landmarks, self.sh)       
		elif self.landmarks.shape != self.sh:
			self.landmarks = self.landmarks.T
		elif self.landmarks.shape == self.sh:
			pass
			#This is the correct shape

		self.n_land = self.landmarks.shape[0]
				
		#Checking index of facets
		if self.facets.min() == 1:
			self.facets = self.facets - 1
			
		return self
	
	def get_normals(self):
		self.normals = get_vertex_normals(self.landmarks, self.facets)
		return self
	
	def get_normal_displacement(self, faceshape):
		#Import libraries
		import numpy as np
		
		diff = self.landmarks - faceshape.landmarks
		self.get_normals()
		normal_displacement = np.squeeze(np.array(np.sum(diff * self.normals, axis = 1)))
		return normal_displacement
	
	def get_euc_dist(self, faceshape):
		'''
		Get the Euclidean distances between to set of 3D landmarks
		Usage
		Input:
			- landmark1 and landmark2: set of 3D landmarks in one single row, in X,Y,Z consecutive fashion, to compute the distances from
		Output:
			- euc_dist: euclidean distance between corresponding set of 3D landmarks
		'''
		#Importing libraries
		import numpy as np
		
		#culating euclidean distance
		euc_dist = np.array(np.sqrt(np.sum(np.power(self.landmarks - faceshape.landmarks, 2), axis=1))).flatten()
		return euc_dist
		
	def take_screenshot(self, profile = 0, colormap = None, colortype = 'GnBu'):
		'''
		Take screenshot of rendered 3D face
		Usage
		Input:
			- landmarks: set of 3D landmarks in one single row, in X,Y,Z consecutive fashion
			- facets: the facets that describe the surface from the 3D landmarks
			- profile: whether to take the screenshot from frontal, midway, or profile (0/1/2)
			- colormap: the colormap to be applied to the 3D surface to the scalar argument
			- colortype: the colormap type
		Output:
			- screenshot: np array of 2D picture
		'''
		#Importing required libraries
		import numpy as np
		from mayavi import mlab

		#Setting offscreen rendering, figure, and X,Y,Z objects
		mlab.options.offscreen = True #True or False
		myfig = mlab.figure(bgcolor=(1, 1, 1), size=(1500, 1500)) 
		X = np.array(self.landmarks[:,0]).flatten()
		Y = np.array(self.landmarks[:,1]).flatten()
		Z = np.array(self.landmarks[:,2]).flatten()
		
		#Generating the mesh, depending on whether colormap is supplied
		if colormap is None:
			mesh = mlab.triangular_mesh(X,Y,Z, self.facets, representation='surface', figure = myfig, color=(0.7,0.7,0.7) )
		elif colormap is not None:
			mesh = mlab.triangular_mesh(X,Y,Z, self.facets, representation='surface', figure = myfig, scalars=colormap, colormap=colortype )
	
		#Setting mesh properties
		mesh.actor.property.backface_culling = True
		mesh.scene.anti_aliasing_frames = 20
		mesh.scene.camera.compute_view_plane_normal()
	
		#Setting camera
		if profile == 0:
			mesh.scene.camera.position = [-0.005, -0.037, 0.067]
			mesh.scene.camera.view_up  = [-0.05, 1, 0.5]
		elif profile == 1:
			mesh.scene.camera.position = [-0.06, -0.03, 0.04]
			mesh.scene.camera.view_up  = [-0.02, 0.87, 0.48]
		elif profile == 2:
			mesh.scene.camera.position = [-0.082, -0.014, 0.0033]
			mesh.scene.camera.view_up  = [-0.12, 0.81, 0.58]
		
		mesh.scene.camera.clipping_range = [0.03, 0.12]
		#Taking screenshot
		screenshot = mlab.screenshot(myfig, mode='rgba', antialiased=True)
		#Set the white background with an alpha of 0, that is, transparent
		screenshot[ np.all(screenshot == 1, axis=2) ] = 0 
		mlab.close()
		return screenshot

	def view_3d(self):
		'''
		Visualize the FaceShape object in 3D using mayavi
		'''

		#Importing required libraries
		import numpy as np
		from mayavi import mlab

		#Setting offscreen rendering, figure, and X,Y,Z objects
		myfig = mlab.figure(bgcolor=(1, 1, 1), size=(1500, 1500))
		X = np.array(self.landmarks[:,0]).flatten()
		Y = np.array(self.landmarks[:,1]).flatten()
		Z = np.array(self.landmarks[:,2]).flatten()

		#Generating the mesh
		mesh = mlab.triangular_mesh(X,Y,Z, self.facets, representation='surface', figure = myfig, color= self.color )

		#Setting mesh properties
		mesh.actor.property.backface_culling = False
		mesh.scene.anti_aliasing_frames = 20
		mesh.scene.camera.compute_view_plane_normal()
		mesh.scene.camera.clipping_range = [0.03, 0.12]
		mesh.scene.camera.position = [-0.005, -0.037, 0.067]
		mesh.scene.camera.view_up  = [-0.05, 1, 0.5]

		mlab.show()
	
def get_vertex_normals(landmarks, facets):
	'''
	Compute the normals at each vertex
	INPUTS:
	 - landmarks: a n(vertices) x 3 array of vertex coordinates or a 3 x n(vertices) x p(observations) array of corresponding vertex co-ordinates
	 - facets:    a n(facets) x 3 connectivity matrix defining the connectivity of the mesh
	OUTPUTS:
	 - Vnormals - a 3 x n(vertices) x n_observations array of vertex normals if if verts_array was three-dimensional, otheriwise it will be a 3 x n(vertices) array of vertex normals
	''' 
	
	import numpy as np
	
	landmarks = np.array(landmarks)
	facets    = np.array(facets)
	
	normals_array = np.zeros_like(landmarks)
	
	N = get_facet_normals(landmarks, facets)
	
	# get edges vectors
	e0 = landmarks[facets[:,2]]-landmarks[facets[:,1]]
	e1 = landmarks[facets[:,0]]-landmarks[facets[:,2]]
	e2 = landmarks[facets[:,1]]-landmarks[facets[:,0]]
	
	# edge length
	de0 = np.linalg.norm(e0,axis=1)
	de1 = np.linalg.norm(e1,axis=1)
	de2 = np.linalg.norm(e2,axis=1)
	
	#face area
	Af = 0.5*np.apply_along_axis(np.linalg.norm,1,np.cross(e0,e1))
	
	Vnormals = np.zeros(landmarks.shape)
	
	for i in range(facets.shape[0]):
		#weight according to area and edge length
		wfv0 = Af[i]/(np.power(de1[i],2)*np.power(de2[i],2))
		wfv1 = Af[i]/(np.power(de0[i],2)*np.power(de2[i],2))
		wfv2 = Af[i]/(np.power(de1[i],2)*np.power(de0[i],2))
	
	
		Vnormals[facets[i,0]] = Vnormals[facets[i,0]] + wfv0*N[i,:]
		Vnormals[facets[i,1]] = Vnormals[facets[i,1]] + wfv1*N[i,:]
		Vnormals[facets[i,2]] = Vnormals[facets[i,2]] + wfv2*N[i,:]
	
	Vnormals = Vnormals/np.atleast_2d(np.linalg.norm(Vnormals,axis=1)).T
	return np.squeeze(Vnormals)

def get_facet_normals(landmarks, facets):
	
	import numpy as np
	
	p0 = landmarks[facets[:,0],:]
	p1 = landmarks[facets[:,1],:]
	p2 = landmarks[facets[:,2],:]
	
	# for each face get two edges
	e0 = p2-p1
	e1 = p0-p2
	
	Fnormals = np.cross(e0,e1)
	
	# make unit length
	Fnormals = (Fnormals.T / np.atleast_2d(np.linalg.norm(Fnormals,axis=1))).T
	return Fnormals
