==================================
Value-View Interface Documentation
==================================

************
Introduction
************

The Value-View Interface (VVI) pattern is a tool to aid in porting complex class structures and design patterns that are not traditionally conducive to execution on the GPU. It can also be used as a stepping stone for performing incremental integration of a large codebase to the GPU. Allowing a code to port it's existing CPU friendly patterns to the GPU without a major initial refactor.

VVI's philosophy is based largely upon ideas from CHAI. VVI builds upon from CHAI's ability to perform implicit data migration between the CPU and GPU for CHAI containers (such as ``chai::ManagedArray``\ ) by extending that capability to user defined types. 

.. note::

   CHAI is designed to integrate as a plugin to RAJA codes, therefore the use of RAJA and CHAI are **essential** for the VVI pattern. 


``chai::ManagedSharedPtr`` (developed specifically for VVI) is the backbone of the VVI Pattern. It is a ``std::shared_ptr`` implementation that supports implicit copies of the owned object between the CPU and GPU as well as dynamic dispatch with base types on the GPU. 

.. note::

   We also provided CHAI with ``chai::ManagedVector``\ , a ``std::vector`` like CHAI class. This was not necessary for VVI but was instrumental in the port for Spheral to the GPU and will be seen later in this document.


*******************
The VALUE Interface
*******************

The Value-View interface enables a class to switch between value (\ *VALUE*\ ) and reference (\ *VIEW*\ )  semantics. A *VALUE* interface allows use of objects in a way that is indicative of normal C++. The semantics of a *VIEW* interface allow us to use the same instance of an object with reference semantics and provides the ability to be implicitly copyable between the HOST and DEVICE.

******************
The VIEW Interface
******************

There are two paths for defining a *VIEW* interface. Each has it's tradeoffs and should be selected based on situation.

Reference Syntax Interface
==========================

The Reference Interface lets a View act like a reference, calling functions of a class similarly to how you would with any reference of a type. This interface is the more verbose and allows users to specify exactly which calls can and cannot be made by the View interface. This is a more type safe implementation. Compile time errors will be thrown if a function call that is only usable in the Value Interface is called from a View Interface. However, if your code targeted for offload to the GPU often references the Implementation class by pointer, this could require many tedious changes and the Pointer Syntax Interface might be a better first step.

Pointer Syntax Interface
========================

The pointer interface allows a View object to act like a pointer. The interface required to define this implementation can be considerably smaller for many classes. It does not provide type safety. Therefore, you *could* call a function that should only be used by a Value object from a View object (you should still get a ``__host__ __device__`` warning if this happens in a GPU context). It is up to the developer to understand what they can and cannot call from a view like object. This is the "quick and dirty" conversion but is really useful for rapid implementation of VVI in certain cases.

*****************
Use Cases for VVI
*****************

VVI was designed to help overcome a variety of complex obstacles in regards to GPU porting, such as:


* Classes with dynamically resizable container members and metadata members associated with the class that need to be consistent across execution contexts.
* Dynamic container-like classes of pointers to other dynamic container-like classes that need to be accessed with double index notation.
* Polymorphic objects with complex implementations. Where they are accessed through base class pointers with dynamic dispatch.
* CRTP compile time polymorphism design patterns across a variety of core classes.

VVI should be avoided when possible. Simple classes with trivial members should be converted to be GPU capable directly. VVI is most useful when a copy of an object could result in a considerable data reallocation or when the class structure is not conducive to traditional GPU use.

***********
Definitions
***********

Implementation Class (\ ``ImplType``\ )
=======================================

The Implementation Class is a class targeted for use with the VVI pattern. This is usually a pre-existing class that a developer wants to use on the GPU. VVI is required when it's current members or class hierarchy / structure are not conventionally suited for device side use.

Value Class (\ ``ValueType``\ )
===============================

An interface to the underlying implmentation class. This class is named the same as the ``ImplType``\ , and thus should require no changes to the current user code to use this type. The ``ValueType`` uses traditional C++ value copy semantics by default. However the semantic interface of how the type should be used should directly match the semantic conditions of the original ``ImplType``.

View Class (\ ``ViewType``\ )
=============================

The View interface is what allows VVI to work seamlessly with ``RAJA`` execution contexts. A ``ViewType`` will have reference semantic like behavior. Meaning that copy construction of a ``ViewType`` does not invoke an allocation of the underlying class data, it instead references the underlying ``ValueType`` instance that it is associated with.

**************************************
Preparing Implementation Class for VVI
**************************************

Private Members
===============


* All data members must be private. They may be accessed through ``get`` and ``set`` methods. This should be done first as this may require some changes through the codebase.

.. important::

   When using VVI the ``ImplType`` members cannot be accessed directly from *VALUE* interfaces so we must expose public members by forwarding a function that returns a reference to them. If we don't make this change, usage of the ``ImplType`` will be different when VVI is enabled vs disabled.


.. code-block:: c++

   class foo {
   public:
     int m_var = 1;
     std::vector<int> m_vec;
   };

   ...

   foo f;
   f.m_var = 5;

.. code-block:: c++

   class foo {
     int m_var = 1;
     std::vector<int> m_vec;
   public:
     int& var() {return m_var}
   };

   ...

   foo f;
   f.var() = 5

std::vector members
===================


* ``std::vector``\ members should be converted to ``vvi::vector``. ``vvi::vector`` is an alias for ``chai::ManagedVector`` that will revert to ``std::vector`` when VVI is disabled.

.. code-block:: c++

   class foo {
     int m_var = 1;
     vvi::vector<int> m_vec;
   public:
     int& var() {return m_var}
   };


* A ``deepCopy()`` method must be provided that calls ``deepCopy()`` on each ``vvi``\ member. When VVI is enabled the ``vvi`` types use reference semantics. ``deepCopy()`` allows the *VALUE* interface to copy the object correctly.

DeepCopy & Compare
  In order to correctly perform Copy and Equivalence operations from the *VALUE* interface, special ``deepCopy`` and ``compare`` method needs to be provided from the ``ImplType``. Helper macros have been developed to aid in declaring these in your classes.

DeepCopy
  ``VVI_IMPL_DEEPCOPY(impl_t, <vvi member names>)`` This deep copy macro first performs the ``ImplType`` Copy Ctor, It will then manually perform a deep-copy on any other ``vvi`` type that exists as a member of the ``ImplType`` defined in the list.

Compare
  ``VVI_IMPL_COMPARE(impl_t, <all class members>)`` This macro will perform a value based comparison on all member names given (you probably want to list all members...).

.. code-block:: c++

   class foo {
     int m_var = 1;
     vvi::vector<int> m_vec;

   public:
     int& var() {return m_var}

     VVI_IMPL_DEEPCOPY(foo, m_vec)
     VVI_IMPL_COMPARE(foo, m_vec, m_var)
   };

chai::CHAICopyable
==================


* If you intend to use the class as the element type of a chai container then a **default constructor** must be provided.
* A class containing a ``vvi`` type member above needs to inherit from ``chai::CHAICopyable``. This allows CHAI to recursively copy the appropriate data .

.. code-block:: c++

   class foo : public chai::CHAICopyable {
     int m_var = 1;
     vvi::vector<int> m_vec;

   public:
     int& var() {return m_var}

     VVI_IMPL_DEEPCOPY(foo, m_vec)
     VVI_IMPL_COMPARE(foo, m_vec, m_var)
   };

   ...

   chai::ManagedArray<foo> my_arr;

Classes w/ virtual methods
==========================


* If the class has a virtual method such that a virtual table would be present for it or its derived classes the class needs to inherit from ``chai::Poly``. 

.. code-block:: c++

   class base { 
     virtual void fooBar() = 0;
   };

   class derived : public base {
     virtual void fooBar() override {...}
   };

.. code-block:: c++

   class base : chai::Poly{ 
     virtual void fooBar() = 0;
   };

   class derived : public base {
     virtual void fooBar() override {...}
   };

.. note::

   Inheriting from ``chai::CHAIPoly`` includes the behavior of ``chai::CHAICopyable``.


*****************************************
Value View Interface Declaration Examples
*****************************************

Basic Class Interface
=====================

Below we demonstrate how to implement the Value-View Interface (VVI) pattern on a very basic class.

We need to wrap our implementation class between ``VVI_IMPL_BEGIN`` and ``VVI_IMPL_END``. These macros encapsulate the code within a namespace ``vvimpl``. This is to allow us to name the future *VALUE* interface class with the same typename as the original implementation.

.. code-block:: c++

   VVI_IMPL_BEGIN
   class foo
   {
     void doSomething() {};
   }
   VVI_IMPL_END

We need to forward declare *VALUE* class names. VVI details should be guarded with a ``VVI_ENABLED`` preprocessor check.

.. code-block:: c++

   #ifdef VVI_ENABLED
   class foo;

Basic VIEW Interface
--------------------

*VIEW* interfaces inherit default constructor, copy constructor, and assignment operator behavior provided by the underlying ``ViewInterface`` type that all *VIEW* interface classes inherit from. 

Most of the Value and View Interface class structure is repeatable and can be rather convoluted to write by hand. ``(PTR/REF)_(VALUE/VIEW)_METACLASS_DECL`` macros are provided to allow users to declare the Value and View interfaces type. 

We use helper macros with the naming convention ``<name>__`` to aid in writing interfaces. Below we pass the type names of ``(<value_type>)``\ , ``(<view_type>)``\ , ``(<impl_type>)`` to each ``_METACLASS_DECL`` macro. The final argument ``(code)`` is used to forward the class definitions further below.

.. note::

   **ALL** arguments to the ``_METACLASS_DECL`` macros **MUST** be wrapped in ``()``.


.. code-block:: c++

   #define fooView__(code) PTR_VIEW_METACLASS_DECL( \
     (foo), (fooView), (vvimpl::foo), (code) ) 

   class fooView__();

.. note::

   We are using the ``PTR_`` implementation of the VVI Pattern in this case. This means all members of the implementation class are available to view objects through ``operator->()``. This leads to many definitions of *VIEW* interfaces to be empty. 

.. tip::

    For empty *VIEW* interfaces, using default behavior as above, we can consolidate the class declaration and definition into one line with ``PTR_VIEW_METACLASS_DEFAULT``

   .. code-block:: c++

       class PTR_VIEW_METACLASS_DEFAULT( (foo), (fooView), (vvimpl::foo) )


Basic Value Interface
---------------------

*VALUE* interface default behavior is inherited from the ``vvi::detail::ValueInterface`` class. If not defined a Default Constructor will be used for the *VALUE* interface. Copy constructor, Assignment and Equivalence behavior is baked in, assuming the appropriate ``DEEPCOPY`` and ``COMPARE`` interfaces have been defined in the ``ImplType``.

VVI Supplies ``VVI_IMPL_INST`` as a way to access the underlying instantiation of the implementation class and should be used to forward functions to the interface as seen below.

.. code-block:: c++

   #define foo__(code) PTR_VALUE_METACLASS_DECL( \
     (foo), (fooView), (vvimpl::foo), (code) )

   class foo__(
   public:
     void doSomething() { VVI_IMPL_INST().doSomething(); }
   );
   #endif

Class w/ Custom Constructor & Vector Member
===========================================

Here we demonstrate using VVI with a class containing a ``vvi::vector`` member and a non default constructor. We have already made the necessary changes for VVI preparation. We will also assume that this class may be used as an element type to another CHAI container, therefore it must inherit from ``chai::CHAICopyable``.

.. code-block:: c++

   VVI_IMPL_BEGIN
   class foo : public chai::CHAICopyable
   {
   public:
     using vec_t = vvi::vector<int>;

     foo(size_t sz) : m_vec(sz) {}

     vec_t& vec() {return m_vec;}
     void doSomething() {};

     VVI_IMPL_DEEPCOPY(foo, m_vec)
   private:
     vec_t m_vec;
   }
   VVI_IMPL_END

Type Alias View Interface
-------------------------

Declaring the interface class names and the *VIEW* interface are the same as above. However we would like to expose the type information for ``vec_t``. Defining the type alias in the *VIEW* interface makes it publicly visible to both *VIEW* and *VALUE* interfaces.

.. code-block:: c++

   #ifdef VVI_ENABLED
   class foo;

   #define fooView__(code) PTR_VIEW_METACLASS_DECL( \
     (foo), (fooView), (vvimpl::foo), (code) ) 

   class fooView__(
   public:
     using vec_t = ImplType::vec_t;
   );

Custom Ctor Value Interface
---------------------------

When defining a non-default constructor for *VALUE* interfaces we need to forward the arguments with ``VVI_VALUE_CTOR_ARGS()`` as below. The argument list needs to be passed within ``()``.

.. code-block:: c++

   #define foo__(code) PTR_VALUE_METACLASS_DECL( \
     (foo), (fooView), (vvimpl::foo), (code) )

   class foo__(
   public:
     foo(size_t sz) : VVI_VALUE_CTOR_ARGS( (sz) ) {}

``vvi::vector`` is an alias to a CHAI container type when VVI is enabled. CHAI containers use reference semantics and need to have their allocations released manually. The containers lifetime should be tied to the lifetime of ``foo``. Therefore, we need to manually free the underlying ``m_vec`` member in the destructor of the value class.

.. code-block:: c++

     ~foo() { VVI_IMPL_INST().m_vec.free(); }

     void doSomething() { VVI_IMPL_INST().doSomething(); }
     vvi::vector<int>& vec() { return VVI_IMPL_INST().vec(); }
   );
   #endif

Template Class Interface
========================

Templated classes can be used with the VVI Pattern.

.. code-block:: c++

   VVI_IMPL_BEGIN
   template<typename T>
   class foo
   {
     T doSomething() {};
   }
   VVI_IMPL_END

Remember to template the forward declaration.

.. code-block:: c++

   #ifdef VVI_ENABLED
   template<typename T>
   class foo;

Template Interface Declarations
-------------------------------

When defining the ``_METACLASS_DECL`` macros we need to template the exclusive types to the interface we are representing:


* For the ``VIEW`` macro: template **only** the ``value`` and ``impl`` type-names.
* For the ``VALUE`` macro : Template **only** the ``view`` and ``impl`` type-names.

.. code-block:: c++

   #define fooView__(code) PTR_VIEW_METACLASS_DECL( \
                           (foo<T>), (fooView), (vvimpl::foo<T>), (code) ) 

   #define foo__(code) PTR_VALUE_METACLASS_DECL( \
                       (foo), (fooView<T>), (vvimpl::foo<T>), (code) )

In the interface declarations we need to ensure we template the interfaces. 

.. code-block:: c++

   template<typename T>
   class fooView__();

   template<typename T>
   class foo__(
   public:
     T doSomething() { return VVI_IMPL_INST().doSomething(); }
   );

   #endif

Abstract & Polymorphic Class Interface
======================================

One of the strongest features of the VVI Pattern is the ability to use abstract class interfaces on the GPU. Below we are converting a *very* basic example of a polymorphic class for use in VVI.

Implmentation Classes w/ Virtual Functions
------------------------------------------

Base classes with virtual interfaces need to inherit from ``vvi::poly``. This allows VVI to perform a specialized copies between the host and device such that the vitual table is preserved across execution spaces.

.. note::

   The implementation classes **MUST** provide default constructors in order to be used with VVI.


.. code-block:: c++

   VVI_IMPL_BEGIN
   class base : public chai::CHAIPoly
   {
     virtual int doSomething() = 0;
   };

   class derived : public base
   {
     derived() : base () {}
     virtual int doSomething() override {}
   };
   VVI_IMPL_END

We need to forward declare for both ``base`` and ``derived``.

.. code-block:: c++

   #ifdef VVI_ENABLED
   class base;
   class derived;

.. code-block:: c++

   #define baseView__(code) PTR_VIEW_METACLASS_DECL( \
     (base), (baseView), (vvimpl::base), (code) ) 
   #define base__(code) PTR_VALUE_METACLASS_DECL( \
     (base), (baseView), (vvimpl::base), (code) )

   #define derivedView__(code) PTR_VIEW_METACLASS_DECL( \
     (derived), (derivedView), (vvimpl::derived), (code) ) 
   #define derived__(code) PTR_VALUE_METACLASS_DECL( \
     (derived), (derivedView), (vvimpl::derived), (code) )

Base Class Interface
--------------------

``base`` is an abstract interface class. It should **NEVER** be constructable. However, we may want to declare the class to query type name information. We use ``VVI_DELETED_INTERFACE()`` to ensure a ``base`` *VALUE* type is not constructed.

.. code-block:: c++

   class baseView__();

   class base__(
     VVI_DELETED_INTERFACE(base) 
   );

   #endif

.. tip::

   There exist shorthand macros for interfaces, the set of which are detailed in sections below.

   .. code-block:: c++

      class PTR_VIEW_METACLASS_DEFAULT( (base), (baseView), (vvimpl::base) )
      class PTR_VALUE_METACLASS_DELETED( (base), (baseView), (vvimpl::base) )


Derived Class Interface
-----------------------

The ``derived`` interface can be defined similarly to other *VALUE* and *VIEW* interfaces. Here we add ``VVI_UPCAST_CONVERSION_OP()`` to allow us to upcast a derived type to the base. This allows conversion from ``derivedView -> baseView`` as well as conversions of ``derived -> baseView``.

.. code-block:: c++

   class derivedView__(
     VVI_UPCAST_CONVERSION_OP(baseView)
   );

   class derived__(
   public:
     T doSomething() { return VVI_IMPL_INST().doSomething(); }
   );

   #endif

*****************
Macro Definitions
*****************

VALUE Interface Definition Macros
=================================

* ``VVI_VALUE_CTOR_ARGS((args))`` : Forward arguments of a non-default constructor for a *VALUE* interface.
* ``VVI_VALUE_DEF_CTOR(value_t)`` : Define default constructor for a *VALUE* interface.


Interface Behavior Macros
=========================

* ``VVI_DELETED_INTERFACE(type)`` : Delete the constructor, copy constructor and assignment operator.
* ``VVI_UPCAST_CONVERSION_OP(parent_t)`` : Define an implicit conversion operator to the defined parent type.


Class Declaration Macros
========================

* ``PTR_VALUE_METACLASS_DECL((value_t), (view_t), (impl_t), (code))``
* ``REF_VALUE_METACLASS_DECL((value_t), (view_t), (impl_t), (code))``
* ``PTR_VIEW_METACLASS_DECL((value_t), (view_t), (impl_t), (code))``
* ``REF_VIEW_METACLASS_DECL((value_t), (view_t), (impl_t), (code))``

* ``PTR_VIEW_METACLASS_DEFAULT((value_t), (view_t), (impl_t))``
* ``PTR_VALUE_METACLASS_DEFAULT((value_t), (view_t), (impl_t))``

* ``PTR_VALUE_METACLASS_DELETED((value_t), (view_t), (impl_t))``
