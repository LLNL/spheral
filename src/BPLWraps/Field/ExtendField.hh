#ifdef __GCCXML__

//------------------------------------------------------------------------------
// Index into a Field.
//------------------------------------------------------------------------------
Scalar getItemScalarField(ScalarField& self, int index);
Vector& getItemVectorField(VectorField& self, int index);
Vector3d& getItemVector3dField(Vector3dField& self, int index);
Tensor& getItemTensorField(TensorField& self, int index);
SymTensor& getItemSymTensorField(SymTensorField& self, int index);
TRTensor& getItemThirdRankTensorField(ThirdRankTensorField& self, int index);
int getItemIntField(IntField& self, int index);
ULL getItemULLField(ULLField& self, int index);
std::vector<double>& getItemVectorDoubleField(VectorDoubleField& self, int index);

//------------------------------------------------------------------------------
// Set a Field value by index.
//------------------------------------------------------------------------------
void setItemScalarField(ScalarField& self, int index, float val);
void setItemVectorField(VectorField& self, int index, Vector& val);
void setItemVector3dField(Vector3dField& self, int index, Vector3d& val);
void setItemTensorField(TensorField& self, int index, Tensor& val);
void setItemSymTensorField(SymTensorField& self, int index, SymTensor& val);
void setItemThirdRankTensorField(ThirdRankTensorField& self, int index, TRTensor& val);
void setItemIntField(IntField& self, int index, int val);
void setItemULLField(ULLField& self, int index, ULL val);
void setItemVectorDoubleField(VectorDoubleField& self, int index, std::vector<double>& val);

//------------------------------------------------------------------------------
// Get slices.
//------------------------------------------------------------------------------
boost::python::list getSliceScalarField(ScalarField& self, int start, int stop);
boost::python::list getSliceVectorField(VectorField& self, int start, int stop);
boost::python::list getSliceVector3dField(Vector3dField& self, int start, int stop);
boost::python::list getSliceTensorField(TensorField& self, int start, int stop);
boost::python::list getSliceSymTensorField(SymTensorField& self, int start, int stop);
boost::python::list getSliceThirdRankTensorField(ThirdRankTensorField& self, int start, int stop);
boost::python::list getSliceIntField(IntField& self, int start, int stop);
boost::python::list getSliceULLField(ULLField& self, int start, int stop);
boost::python::list getSliceVectorDoubleField(VectorDoubleField& self, int start, int stop);

//------------------------------------------------------------------------------
// Set slices.
//------------------------------------------------------------------------------
void setSliceScalarField(ScalarField& self, int start, int stop, const boost::python::list& values);
void setSliceVectorField(VectorField& self,int start, int stop, const boost::python::list& values);
void setSliceVector3dField(Vector3dField& self, int start, int stop, const boost::python::list& values);
void setSliceTensorField(TensorField& self, int start, int stop, const boost::python::list& values);
void setSliceSymTensorField(SymTensorField& self, int start, int stop, const boost::python::list& values);
void setSliceThirdRankTensorField(ThirdRankTensorField& self, int start, int stop, const boost::python::list& values);
void setSliceIntField(IntField& self, int start, int stop, const boost::python::list& values);
void setSliceULLField(ULLField& self, int start, int stop, const boost::python::list& values);
void setSliceVectorDoubleField(VectorDoubleField& self, int start, int stop, const boost::python::list& values);

//------------------------------------------------------------------------------
// Provide a list interface for the (all, internal, ghost) values.
//------------------------------------------------------------------------------
boost::python::list allScalarValues(ScalarField* self);
boost::python::list internalScalarValues(ScalarField* self);
boost::python::list ghostScalarValues(ScalarField* self);

boost::python::list allVectorValues(VectorField* self);
boost::python::list internalVectorValues(VectorField* self);
boost::python::list ghostVectorValues(VectorField* self);

boost::python::list allVector3dValues(Vector3dField* self);
boost::python::list internalVector3dValues(Vector3dField* self);
boost::python::list ghostVector3dValues(Vector3dField* self);

boost::python::list allTensorValues(TensorField* self);
boost::python::list internalTensorValues(TensorField* self);
boost::python::list ghostTensorValues(TensorField* self);

boost::python::list allSymTensorValues(SymTensorField* self);
boost::python::list internalSymTensorValues(SymTensorField* self);
boost::python::list ghostSymTensorValues(SymTensorField* self);

boost::python::list allThirdRankTensorValues(ThirdRankTensorField* self);
boost::python::list internalThirdRankTensorValues(ThirdRankTensorField* self);
boost::python::list ghostThirdRankTensorValues(ThirdRankTensorField* self);

boost::python::list allIntValues(IntField* self);
boost::python::list internalIntValues(IntField* self);
boost::python::list ghostIntValues(IntField* self);

boost::python::list allULLValues(ULLField* self);
boost::python::list internalULLValues(ULLField* self);
boost::python::list ghostULLValues(ULLField* self);

boost::python::list allVectorDoubleValues(VectorDoubleField* self);
boost::python::list internalVectorDoubleValues(VectorDoubleField* self);
boost::python::list ghostVectorDoubleValues(VectorDoubleField* self);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#else
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// Index into a Field.
//------------------------------------------------------------------------------
inline
Scalar
getItemScalarField(ScalarField& self,
                   int index) {
  return getItemByValue<ScalarField, Scalar>(self, index);
}

inline
Vector& 
getItemVectorField(VectorField& self,
                   int index) {
  return getItemByReference<VectorField, Vector>(self, index);
}

inline
Vector3d& 
getItemVector3dField(Vector3dField& self,
                     int index) {
  return getItemByReference<Vector3dField, Vector3d>(self, index);
}

inline
Tensor& 
getItemTensorField(TensorField& self,
                   int index) {
  return getItemByReference<TensorField, Tensor>(self, index);
}

inline
SymTensor& 
getItemSymTensorField(SymTensorField& self,
                      int index) {
  return getItemByReference<SymTensorField, SymTensor>(self, index);
}

inline
TRTensor& 
getItemThirdRankTensorField(ThirdRankTensorField& self,
                      int index) {
  return getItemByReference<ThirdRankTensorField, TRTensor>(self, index);
}

inline
int
getItemIntField(IntField& self,
                int index) {
  return getItemByValue<IntField, ULL>(self, index);
}

inline
int
getItemULLField(ULLField& self,
                int index) {
  return getItemByValue<ULLField, int>(self, index);
}

inline
std::vector<double>&
getItemVectorDoubleField(VectorDoubleField& self,
                         int index) {
  return getItemByReference<VectorDoubleField, std::vector<double> >(self, index);
}

//------------------------------------------------------------------------------
// Set a Field value by index.
//------------------------------------------------------------------------------
inline
void
setItemScalarField(ScalarField& self,
                   int index,
                   float val) {
  setItemByValue<ScalarField, Scalar>(self, index, val);
}

inline
void
setItemVectorField(VectorField& self,
                   int index,
                   Vector& val) {
  setItemByReference<VectorField, Vector>(self, index, val);
}

inline
void
setItemVector3dField(Vector3dField& self,
                     int index,
                     Vector3d& val) {
  setItemByReference<Vector3dField, Vector3d>(self, index, val);
}

inline
void
setItemTensorField(TensorField& self,
                   int index,
                   Tensor& val) {
  setItemByReference<TensorField, Tensor>(self, index, val);
}

inline
void
setItemSymTensorField(SymTensorField& self,
                      int index,
                      SymTensor& val) {
  setItemByReference<SymTensorField, SymTensor>(self, index, val);
}

inline
void
setItemThirdRankTensorField(ThirdRankTensorField& self,
                      int index,
                      TRTensor& val) {
  setItemByReference<ThirdRankTensorField, TRTensor>(self, index, val);
}

inline
void
setItemIntField(IntField& self,
                int index,
                int val) {
  setItemByValue<IntField, ULL>(self, index, val);
}

inline
void
setItemULLField(ULLField& self,
                int index,
                ULL val) {
  setItemByValue<ULLField, ULL>(self, index, val);
}

inline
void
setItemVectorDoubleField(VectorDoubleField& self,
                         int index,
                         std::vector<double>& val) {
  setItemByReference<VectorDoubleField, std::vector<double> >(self, index, val);
}

//------------------------------------------------------------------------------
// Get slices.
//------------------------------------------------------------------------------
inline
boost::python::list
getSliceScalarField(ScalarField& self,
                    int start,
                    int stop) {
  return getSlice<ScalarField>(self, start, stop);
}

inline
boost::python::list
getSliceVectorField(VectorField& self,
                    int start,
                    int stop) {
  return getSlice<VectorField>(self, start, stop);
}

inline
boost::python::list
getSliceVector3dField(Vector3dField& self,
		      int start,
		      int stop) {
  return getSlice<Vector3dField>(self, start, stop);
}

inline
boost::python::list
getSliceTensorField(TensorField& self,
                    int start,
                    int stop) {
  return getSlice<TensorField>(self, start, stop);
}

inline
boost::python::list
getSliceSymTensorField(SymTensorField& self,
                       int start,
                       int stop) {
  return getSlice<SymTensorField>(self, start, stop);
}

inline
boost::python::list
getSliceThirdRankTensorField(ThirdRankTensorField& self,
                       int start,
                       int stop) {
  return getSlice<ThirdRankTensorField>(self, start, stop);
}

inline
boost::python::list
getSliceIntField(IntField& self,
                 int start,
                 int stop) {
  return getSlice<IntField>(self, start, stop);
}

inline
boost::python::list
getSliceULLField(ULLField& self,
                 int start,
                 int stop) {
  return getSlice<ULLField>(self, start, stop);
}

inline
boost::python::list
getSliceVectorDoubleField(VectorDoubleField& self,
                          int start,
                          int stop) {
  return getSlice<VectorDoubleField>(self, start, stop);
}

//------------------------------------------------------------------------------
// Set slices.
//------------------------------------------------------------------------------
inline
void
setSliceScalarField(ScalarField& self,
                    int start,
                    int stop,
                    const boost::python::list& values) {
  return setSlice<ScalarField, Scalar>(self, start, stop, values);
}

inline
void
setSliceVectorField(VectorField& self,
                    int start,
                    int stop,
                    const boost::python::list& values) {
  return setSlice<VectorField, Vector>(self, start, stop, values);
}

inline
void
setSliceVector3dField(Vector3dField& self,
		      int start,
		      int stop,
		      const boost::python::list& values) {
  return setSlice<Vector3dField, Vector3d>(self, start, stop, values);
}

inline
void
setSliceTensorField(TensorField& self,
                    int start,
                    int stop,
                    const boost::python::list& values) {
  return setSlice<TensorField, Tensor>(self, start, stop, values);
}

inline
void
setSliceSymTensorField(SymTensorField& self,
                       int start,
                       int stop,
                       const boost::python::list& values) {
  return setSlice<SymTensorField, SymTensor>(self, start, stop, values);
}

inline
void
setSliceThirdRankTensorField(ThirdRankTensorField& self,
                       int start,
                       int stop,
                       const boost::python::list& values) {
  return setSlice<ThirdRankTensorField, TRTensor>(self, start, stop, values);
}

inline
void
setSliceIntField(IntField& self,
                 int start,
                 int stop,
                 const boost::python::list& values) {
  return setSlice<IntField, int>(self, start, stop, values);
}

inline
void
setSliceULLField(ULLField& self,
                 int start,
                 int stop,
                 const boost::python::list& values) {
  return setSlice<ULLField, ULL>(self, start, stop, values);
}

inline
void
setSliceVectorDoubleField(VectorDoubleField& self,
                          int start,
                          int stop,
                          const boost::python::list& values) {
  return setSlice<VectorDoubleField, std::vector<double> >(self, start, stop, values);
}

//------------------------------------------------------------------------------
// Provide a list interface for the (all, internal, ghost) values.
//------------------------------------------------------------------------------
// Scalar Field
inline
boost::python::list
allScalarValues(ScalarField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalScalarValues(ScalarField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostScalarValues(ScalarField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// Vector Field
inline
boost::python::list
allVectorValues(VectorField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalVectorValues(VectorField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostVectorValues(VectorField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// Vector3d Field
inline
boost::python::list
allVector3dValues(Vector3dField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalVector3dValues(Vector3dField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostVector3dValues(Vector3dField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// Tensor Field
inline
boost::python::list
allTensorValues(TensorField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalTensorValues(TensorField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostTensorValues(TensorField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// SymTensor Field
inline
boost::python::list
allSymTensorValues(SymTensorField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalSymTensorValues(SymTensorField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostSymTensorValues(SymTensorField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// ThirdRankTensor Field
inline
boost::python::list
allThirdRankTensorValues(ThirdRankTensorField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalThirdRankTensorValues(ThirdRankTensorField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostThirdRankTensorValues(ThirdRankTensorField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// Int Field
inline
boost::python::list
allIntValues(IntField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalIntValues(IntField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostIntValues(IntField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// ULL Field
inline
boost::python::list
allULLValues(ULLField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalULLValues(ULLField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostULLValues(ULLField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

// VectorDouble Field
inline
boost::python::list
allVectorDoubleValues(VectorDoubleField* self) {
  return Spheral::iteratorsAsListByValue(self->begin(), 
                                         self->end());
}

inline
boost::python::list
internalVectorDoubleValues(VectorDoubleField* self) {
  return Spheral::iteratorsAsListByValue(self->internalBegin(), 
                                         self->internalEnd());
}

inline
boost::python::list
ghostVectorDoubleValues(VectorDoubleField* self) {
  return Spheral::iteratorsAsListByValue(self->ghostBegin(), 
                                         self->ghostEnd());
}

#endif
