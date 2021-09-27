from Spheral import *
import numpy as np

# 09/27/2020 -- generalized for databases consisting of Fluid/Solid node lists.
# general wrapper for different sampling objects
class Sample:
    def __init__(self,db,sampleFunc,samplers,dataClipper=None):
        self.db = db
        self.samplers = samplers
        self.sampleFunc = sampleFunc

        if dataClipper:
          self.dataClipper = dataClipper
        else:
          self.dataClipper = self.dummyTrueFunc
    
    def dummyTrueFunc(self,pos): 
        return True

    def sample(self,cycle,t,dt):

        pos = self.db.globalPosition

        # clear out anything from last sampling cycle
        for sampler in self.samplers:
            sampler.clear()

        # loop over node lists, nodes, sample functions
        for j in range(pos.numFields):

            posj = pos[j]
            ni = posj.numInternalElements
            for i in range(ni) :
                if self.dataClipper(posj[i]):
                    for sampler in self.samplers:
                        sampler.sampleNode(posj[i],self.sampleFunc(j,i))
        
        # write the data   
        for sampler in self.samplers:
            sampler.write(cycle,t)



#===================================================================
# Base sample type holding file io and other generic attributes
#===================================================================
class SampleTypeBase:
    def __init__(self,
                 N_Fields, # number of fields being sampled (theres a better way to do this)
                 histDir, # Directory for output files
                 n_samples,
                 dataClipper=None, # function to clip the data 
                 header_label = ' ', # header for output file
                 WT=TableKernel3d(NBSplineKernel3d(7), 100), # SPH smoothing kernel used to sample
                 ):
                 
        self.WT = WT
        self.histDir = histDir
        self.header_label = header_label
        self.kernelFactor = WT.kernelExtent/2.0 # normalize to span 2 neighbor samples

        if dataClipper:
          self.dataClipper = dataClipper
        else:
          def dataClipper_func(pos): 
            return True
          self.dataClipper=dataClipper_func

        # preallocating data output
        self.DATA = np.zeros((n_samples,N_Fields))
        self.wsum = np.zeros((n_samples,N_Fields))
        self.DATA_Tot = np.zeros((n_samples,N_Fields+3))


    def clear(self): # freshen things up for next round of sampling
    
        self.DATA[:,:]=0.0
        self.wsum[:,:]=0.0
        self.DATA_Tot[:,3:]=0.0 # dont clear positions

    def write(self,cycle,t): # dump to gnu files
        
        self.DATA_Tot[:,3:]=mpi.allreduce(self.DATA,mpi.SUM)/mpi.allreduce(self.wsum,mpi.SUM) 
        np.nan_to_num(self.DATA_Tot[:,3:], copy=False) 

        if mpi.rank == 0:
            np.savetxt(self.histDir+("cycle=%010i_time=%gsec.gnu" % (cycle, t)),self.DATA_Tot,delimiter=' ',header=self.header_label)


#===================================================================
# Line Segment in 2d or 3d defined by two points
#===================================================================
class LineSegment(SampleTypeBase):

    def __init__(self,
                 p0, # line segment start coordinates
                 p1, # line segment end coordintes
                 resolution, # resolution of sample points along line
                 N_Fields, # number of fields being sampled (theres a better way to do this)
                 histDir, # Directory for output files
                 dataClipper=None, # function to clip the data 
                 header_label = None, # header for output file
                 WT=TableKernel3d(NBSplineKernel3d(7), 100), # SPH smoothing kernel used to sample
                 ):
        
        
        self.p0 = p0
        self.n_samples = max(int((p1-p0).magnitude()/resolution),2)

        SampleTypeBase.__init__(self,
                 N_Fields, # number of fields being sampled (theres a better way to do this)
                 histDir, # Directory for output files
                 self.n_samples,
                 dataClipper = dataClipper, # function to clip the data 
                 header_label = header_label, # header for output file
                 WT = WT)

        
        # tangent coordinates of segment
        self.tan_vec=(p1-self.p0).unitVector()
        self.S0 = self.p0.dot(self.tan_vec)
        self.S1 = p1.dot(self.tan_vec)
        self.Smax = max(self.S0,self.S1) 
        self.Smin = min(self.S0,self.S1)
        self.segmentLength=abs(self.S1-self.S0)

        sample_x = np.linspace(0.0,1.0,self.n_samples) #normalized sample points
        
        self.res = (self.Smax-self.Smin)/(self.n_samples-1) # resolution of line

        # set the positions of the sample points in the data output file
        self.DATA_Tot[:,0]=(p1.x-self.p0.x)*sample_x.transpose()+self.p0.x
        self.DATA_Tot[:,1]=(p1.y-self.p0.y)*sample_x.transpose()+self.p0.y 
        if p1.z:
            self.DATA_Tot[:,2]=(p1.z-self.p0.z)*sample_x.transpose()+self.p0.z 

    def sampleNode(self,pos,Fields):
       
        pos_offset = pos-self.p0 # set starting point as origin
        si = pos_offset.dot(self.tan_vec) # tangent coordinate
        ni = (pos_offset.cross(self.tan_vec)).magnitude()/self.res # normal coordinate

        # issue w/ if statemennt and min max stuff, needs to be cleaned up
        # test if point is within segment bounds along tangent and within range
        # of the kernel in the normal. Apply clipping test too.
        if (si<self.segmentLength and si>0.0) and ni < 2.0 and self.dataClipper(pos):
            
            index_float = si/self.segmentLength*(self.n_samples-1)
            index_int = int(index_float) - 1
            interp = index_float-index_int + 1
            Vars = np.array(Fields) 
            
            for i in range(4):
                if (index_int+i)>0 and (index_int+i)<self.n_samples:
                    ds = abs(interp - i)
                    W = self.WT.kernelValue(self.kernelFactor*Vector2d(ni,ds).magnitude(),1.0)**2
                    self.DATA[index_int+i,:] += W*Vars 
                    self.wsum[index_int+i,:] += W 



#===================================================================
# Plane in 3d defined by point and 3 vectors
#===================================================================
class Rectangle(SampleTypeBase):

    def __init__(self,
                 p0, # central point (anchor of plane)
                 p1, # establishes direction 1 
                 p2, # establishes direction 2
                 resolution, # resolution of sample points along line
                 N_Fields, # number of fields being sampled (theres a better way to do this)
                 histDir, # Directory for output files
                 dataClipper=None, # function to clip the data 
                 header_label = None, # header for output file
                 WT=TableKernel3d(NBSplineKernel3d(7), 100), # SPH smoothing kernel used to sample
                 ):
        
        
        self.p0 = p0
        self.n = (p1.cross(p2)).unitVector()
        
        self.N_samples_u = max(int((p1-p0).magnitude()/resolution),2)
        self.N_samples_v = max(int((p2-p0).magnitude()/resolution),2)
        
        self.N = self.N_samples_u*self.N_samples_v
        
        SampleTypeBase.__init__(self,
                 N_Fields, # number of fields being sampled (theres a better way to do this)
                 histDir, # Directory for output files
                 self.n_samples,
                 dataClipper = dataClipper, # function to clip the data 
                 header_label = header_label, # header for output file
                 WT = WT)

        
        # in-plane direction 1
        self.tan_vec_u =(p1-self.p0).unitVector()
        self.S0_u = self.p0.dot(self.tan_vec_u)
        self.S1_u = p1.dot(self.tan_vec_u)
        self.Smax_u = max(self.S0_u,self.S1_u) 
        self.Smin_u= min(self.S0_u,self.S1_u)

        # in-plane direction 2
        self.tan_vec_v =(p2-self.p0).unitVector()
        self.S0_v = self.p0.dot(self.tan_vec_v)
        self.S1_v = p2.dot(self.tan_vec_v)
        self.Smax_v = max(self.S0_v,self.S1_v) 
        self.Smin_v= min(self.S0_v,self.S1_v)
        
        # max extent for clipping later
        res_u = (self.Smax_u-self.Smin_u)/(self.N_samples_u-1) # resolution of line
        res_v = (self.Smax_v-self.Smin_v)/(self.N_samples_v-1) # resolution of line
        self.res = max(res_u,res_v)
        
        # set the positions of the sample points in the data output file
        iter=0
        for i in range(self.N_samples_u):
            for j in range(self.N_samples_v):
                
                coordi = float(i)/float(self.N_samples_u-1) # compatibility w/ python 2
                coordj = float(j)/float(self.N_samples_v-1)

                pos_sample = coordi*(p1-self.p0) + jcoord*(p2-self.p0) + self.p0
               
                self.DATA_Tot[iter,0]=pos_sample[0]
                self.DATA_Tot[iter,1]=pos_sample[1]
                self.DATA_Tot[iter,2]=pos_sample[2]

    def sampleNode(self,pos,Fields):

        pos_offset = pos-self.p0 # set starting point as origin
        si_u = pos_offset.dot(self.tan_vec_u) # tangent coordinate
        si_v = pos_offset.dot(self.tan_vec_v) # tangent coordinate
        ni = (pos_offset.dot(self.n)).magnitude()/self.res # normal coordinate

        bounds = (si_u<self.Smax_u and si_u>self.Smin_u) and (si_v<self.Smax_v and si_v>self.Smin_v)
        
        
        if bounds and ni < 2.0 and self.dataClipper(pos):
        
            Vars = np.array(Fields)
            
            # direction 1
            index_float_u = (si_u-self.S0_u)/(self.S1_u-self.S0_u)*(self.N_samples_u-1)
            index_int_u = int(index_float_u) 
            interp_u = index_float_u-index_int_u
            
            # direction 2
            index_float_v = (si_v-self.S0_v)/(self.S1_v-self.S0_v)*(self.N_samples_v-1)
            index_int_v = int(index_float_v) 
            interp_v = index_float_v-index_int_v
            
            # base index
            index_base = self.N_samples_v*(index_int_u-1) + (index_int_v-1)
            ds_u = interp_u+1
            ds_v = interp_v+1
            
            # loop through the 16 sample points that 
            # could potentially interact with the particle
            for i in range(16):
                index = index_base + self.N_samples_v*i/4 + i%4 # index
                if index >= 0  and index < self.N: # make sure were in range
                    ds = np.sqrt((ds_u-i/4)**2+(ds_v-i%4)**2)
                    W = self.WT.kernelValue(self.kernelFactor*Vector2d(ni,ds).magnitude(),1.0)**2
                    self.DATA[index,:] += W*Vars 
                    self.wsum[index,:] += W 
                
        

