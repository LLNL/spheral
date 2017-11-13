//------------------------------------------------------------------------------
// NodeListRegistrar
// 
// A singleton class that maintains a running (sorted) record of NodeList names 
// that are currently in scope.  This is necessary to ensure that there is a
// globally available unique order to iterate over NodeLists (and Fields in 
// FieldLists) that will be uniform across processors.
//
// Created by JMO, Fri Aug  5 10:16:25 PDT 2005
//------------------------------------------------------------------------------
#ifndef __Spheral_NodeListRegistrar__
#define __Spheral_NodeListRegistrar__

#include <string>

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

namespace Spheral {

// Forward decalarations.
namespace NodeSpace {
  template<typename Dimension> class NodeList;
  template<typename Dimension> class FluidNodeList;
}
namespace FieldSpace {
  template<typename Dimension> class FieldBase;
}

template<typename Dimension>
class NodeListRegistrar {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef std::vector<NodeSpace::NodeList<Dimension>*> ContainerType;
  typedef std::vector<NodeSpace::FluidNodeList<Dimension>*> FluidContainerType;

  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;

  typedef typename FluidContainerType::iterator fluid_iterator;
  typedef typename FluidContainerType::const_iterator const_fluid_iterator;

  // Define a nested comparator class for comparing NodeLists by their names.
  class NodeListComparator {
  public:
    int operator()(const NodeSpace::NodeList<Dimension>* nodeListPtr1,
                   const NodeSpace::NodeList<Dimension>* nodeListPtr2) const {
      return nodeListPtr1->name() < nodeListPtr2->name();
    }
    int operator()(const FieldSpace::FieldBase<Dimension>* fieldPtr1,
                   const FieldSpace::FieldBase<Dimension>* fieldPtr2) const {
      return fieldPtr1->nodeListPtr()->name() < fieldPtr2->nodeListPtr()->name();
    }
  };

  // Get the instance.
  static NodeListRegistrar& instance();
//#pragma omp declare target
  static NodeListRegistrar& getInstance();
//#pragma omp end declare target
  // The number of NodeLists currently in play.
  int numNodeLists() const;
  int numFluidNodeLists() const;

  // Iterators over the current NodeLists (in proper order of course).
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  fluid_iterator fluidBegin();
  fluid_iterator fluidEnd();

  const_fluid_iterator fluidBegin() const;
  const_fluid_iterator fluidEnd() const;

  // The (sorted) set of registered NodeList names.
  std::vector<std::string> registeredNames() const;
  std::vector<std::string> registeredFluidNames() const;

  // Internal consistency checking.
  bool valid() const;

  // Determine the proper place in a sequence of Fields that a given Field
  // should be inserted.
  // This is primarily a helper for building FieldLists properly.
  template<typename IteratorType, typename ThingyType>
  IteratorType
  findInsertionPoint(const ThingyType& thingy,
                     const IteratorType begin,
                     const IteratorType end) const;

  // Flag (for use by others) indicating whether we want to run in 
  // strictly domain decomposition independent/reproducing mode (may have
  // a performance impact).
  bool domainDecompositionIndependent() const;
//#pragma omp declare target
  void domainDecompositionIndependent(const bool x);
//#pragma omp end declare target
private:
  //--------------------------- Private Interface---------------------------//
#ifndef __GCCXML__
  // The one and only instance.
//#pragma omp declare target
  static NodeListRegistrar* mInstancePtr;
//#pragma omp end declare target
  // The current set of NodeLists.
  ContainerType mNodeLists;
  FluidContainerType mFluidNodeLists;

  // Flag for the domain independent choice.
//#pragma omp declare target
  bool mDomainDecompIndependent;
//#pragma omp end declare target
  // No public constructors, destructor, or assignment.
  NodeListRegistrar();
  NodeListRegistrar(const NodeListRegistrar&);
  NodeListRegistrar& operator=(const NodeListRegistrar&);
  ~NodeListRegistrar();

  // Grant friendship to NodeLists to register and unregister themselves.
  friend class NodeSpace::NodeList<Dimension>;
  friend class NodeSpace::FluidNodeList<Dimension>;
  void registerNodeList(NodeSpace::NodeList<Dimension>& nodeList);
  void registerNodeList(NodeSpace::FluidNodeList<Dimension>& nodeList);
  void unregisterNodeList(NodeSpace::NodeList<Dimension>& nodeList);
  void unregisterNodeList(NodeSpace::FluidNodeList<Dimension>& nodeList);

  // Helpers for findInsertionPoint that know how to extract a NodeList*
  // from the argument.
  NodeSpace::NodeList<Dimension>* 
  getNodeListPtr(NodeSpace::NodeList<Dimension>* thingy) const {
    return thingy;
  }

  NodeSpace::NodeList<Dimension>* 
  getNodeListPtr(const NodeSpace::NodeList<Dimension>* thingy) const {
    return const_cast<NodeSpace::NodeList<Dimension>*>(thingy);
  }

  NodeSpace::NodeList<Dimension>* 
  getNodeListPtr(NodeSpace::FluidNodeList<Dimension>* thingy) const {
    return (NodeSpace::NodeList<Dimension>*) thingy;
  }

  NodeSpace::NodeList<Dimension>* 
  getNodeListPtr(FieldSpace::FieldBase<Dimension>* thingy) const {
    return const_cast<NodeSpace::NodeList<Dimension>*>(thingy->nodeListPtr());
  }

  NodeSpace::NodeList<Dimension>* 
  getNodeListPtr(const FieldSpace::FieldBase<Dimension>* thingy) const {
    return const_cast<NodeSpace::NodeList<Dimension>*>(thingy->nodeListPtr());
  }
#endif

};
}

#ifndef __GCCXML__
#include "NodeListRegistrarInline.hh"
#endif
#else

// Forward declaration.

namespace Spheral {
  template<typename Dimension>  class NodeListRegistrar;
}
#endif
