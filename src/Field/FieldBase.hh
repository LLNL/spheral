//---------------------------------Spheral++----------------------------------//
// FieldBase -- An abstract base class to provide a generic handle on Fields.
//
// Created by JMO, Tue Sep 14 22:36:19 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldBase_hh__
#define __Spheral_FieldBase_hh__

#include <string>
#include <vector>
#include <memory>

#include "Utilities/ValueViewInterface.hh"

namespace Spheral {

// Forward declarations
template<typename Dimension> class FieldListBase;
template<typename Dimension> class NodeList;

VVI_IMPL_BEGIN

template<typename Dimension>
class FieldBase : chai::CHAIPoly {
public:
  //--------------------------- Public Interface ---------------------------//
  using FieldName = std::string;
  typedef typename Dimension::Scalar Scalar;

protected:

  //--------------------------- Private Interface ---------------------------//
  FieldName mName;
  const NodeList<Dimension>* mNodeListPtr = 0;

  // Disallow the default constructor.
  FieldBase();
  // Constructors.
  FieldBase(FieldName name);
  FieldBase(FieldName name,
            const NodeList<Dimension>& nodeList);

  FieldBase(const FieldBase& fieldBase);
  // Assignment operator.
  virtual FieldBase& operator=(const FieldBase& rhs) = default;


public:
  //--------------------------- Public Interface ---------------------------//
  virtual std::shared_ptr<FieldBase> clone() const = 0;

  // Destructor.
  virtual ~FieldBase();

  // Require descendent fields be able to test equivalence.
  virtual bool operator==(const FieldBase& rhs) const = 0;
  bool operator!=(const FieldBase& rhs) const;

  // Access the name.
  FieldName name() const;
  void name(FieldName name);

  // Provide methods to access and set the NodeList.
  const NodeList<Dimension>& nodeList() const;
  const NodeList<Dimension>* nodeListPtr() const;
  void unregisterNodeList();

  // Methods every field must provide.
  virtual unsigned size() const = 0;
  virtual void Zero() = 0;
  virtual void setNodeList(const NodeList<Dimension>& nodeList) = 0;
  virtual void resizeField(unsigned size) = 0;
  virtual void resizeFieldInternal(unsigned size, unsigned oldFirstGhostNode) = 0;
  virtual void resizeFieldGhost(unsigned size) = 0;
  virtual void deleteElement(int nodeID) = 0;
  virtual void deleteElements(const std::vector<int>& nodeIDs) = 0;
  virtual std::vector<char> packValues(const std::vector<int>& nodeIDs) const = 0;
  virtual void unpackValues(const std::vector<int>& nodeIDs,
                            const std::vector<char>& buffer) = 0;
  virtual void copyElements(const std::vector<int>& fromIndices,
                            const std::vector<int>& toIndices) = 0;
  virtual bool fixedSizeDataType() const = 0;
  virtual int numValsInDataType() const = 0;
  virtual int sizeofDataType() const = 0;
  virtual int computeCommBufferSize(const std::vector<int>& packIndices,
                                    const int sendProc,
                                    const int recvProc) const = 0;

protected:
  //--------------------------- Protected Interface ---------------------------//
  void setFieldBaseNodeList(const NodeList<Dimension>& nodeListPtr);

};

VVI_IMPL_END

#ifdef VVI_ENABLED
//-----------------------------------------------------------------------------
// Interface to support porting to the GPU.
//-----------------------------------------------------------------------------

// We need to forward declare value classes for view interface definitions.
template<typename Dimension>
class FieldBase;

// Define Metaclass macros for Value/View relationships
template<typename Dimension>
class PTR_VIEW_METACLASS_DEFAULT((FieldBase<Dimension>), (FieldBaseView), (vvimpl::FieldBase<Dimension>))

template<typename Dimension>
class PTR_VALUE_METACLASS_DELETED((FieldBase), (FieldBaseView<Dimension>), (vvimpl::FieldBase<Dimension>))

#endif // !defined(SPHERAL_ENABLE_VVI)

} // namespace Spheral

#include "FieldBaseInline.hh"

#endif
