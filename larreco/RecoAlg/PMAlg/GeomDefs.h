/**
 * @file   GeomDefs.h
 * @brief  Definition of some geometry-related data structures
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 8, 2016
 * 
 * This is a pure header: no implementation file is required.
 */

#ifndef LARRECO_RECOALG_PMALG_GEOMDEFS_H
#define LARRECO_RECOALG_PMALG_GEOMDEFS_H

// ROOT
#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/PositionVector2D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/DisplacementVector2D.h"
#include "Math/GenVector/DisplacementVector3D.h"

// C/C++ standard libraries
#include <type_traits>
#include <iterator> // std::begin, std::end
#include <utility> // std::forward()
#include <cassert> // assert()


namespace pma {
  
  /// Type of real number used in the computations (double precision)
  using Real = double;
  
  /// coordinate system used for our 3D objects
  using CoordSystem3D = ROOT::Math::Cartesian3D<Real>;
  
  /// coordinate system used for our 2D objects
  using CoordSystem2D = ROOT::Math::Cartesian2D<Real>;
  
  /// A point in 2D space
  using Point2D_t = ROOT::Math::PositionVector2D<CoordSystem2D>;
  
  /// A point in 3D space
  using Point3D_t = ROOT::Math::PositionVector3D<CoordSystem3D>;
  
  
  /// A vector in 2D space
  using Vector2D_t = ROOT::Math::DisplacementVector2D<CoordSystem2D>;
  
  /// A vector in 3D space
  using Vector3D_t = ROOT::Math::DisplacementVector3D<CoordSystem3D>;
  
  
  template <typename Vect>
  class PointAverage;
  
  
} // namespace pma


//------------------------------------------------------------------------------
//--- template implementation
//---

namespace pma {
  namespace details {
    
    //--------------------------------------------------------------------------
    template <typename Vect, unsigned int N, typename = void>
    struct MinVectorDim: public std::false_type
      { static_assert(N <= 4, "MinVectorDim supports at most dimension 3."); };

    template <typename Vect, typename Dummy>
    struct MinVectorDim<Vect, 0, Dummy>: public std::true_type {};

    template <typename Vect>
    struct MinVectorDim<
      Vect, 1U,
      std::enable_if_t
        <std::is_arithmetic<std::decay_t<decltype(Vect().X())>>::value, void>
      >
      : public std::true_type {};

    template <typename Vect>
    struct MinVectorDim<
      Vect, 2U,
      std::enable_if_t
        <std::is_arithmetic<std::decay_t<decltype(Vect().Y())>>::value, void>
      >
      : public MinVectorDim<Vect, 1U> {};

    template <typename Vect>
    struct MinVectorDim<
      Vect, 3U,
      std::enable_if_t
        <std::is_arithmetic<std::decay_t<decltype(Vect().Z())>>::value, void>
      >
      : public MinVectorDim<Vect, 2U> {};

    template <typename Vect>
    struct MinVectorDim<
      Vect, 4U,
      std::enable_if_t
        <std::is_arithmetic<std::decay_t<decltype(Vect().T())>>::value, void>
      >
      : public MinVectorDim<Vect, 3U> {};


    template <typename Vect, unsigned int N>
    struct VectorDimNotAbove
      : public std::integral_constant<bool, !MinVectorDim<Vect, N+1>::value>
    {
      static constexpr unsigned int Dim = N;
    };


    template <typename Vect, unsigned int Min = 0, typename = void>
    struct MaxVectorDim {
      static constexpr unsigned int value = Min;
    };

    template <typename Vect, unsigned int Min>
    struct MaxVectorDim
      <Vect, Min, std::enable_if_t<MinVectorDim<Vect, (Min+1)>::value, void>>
      : public MaxVectorDim<Vect, Min + 1>
    {};

    template <typename Vect>
    using VectorDim = MaxVectorDim<Vect>;
    
    
    //--------------------------------------------------------------------------
    template <typename Vect, unsigned int Dim>
    struct VectorCoordStruct;
    
    template <typename Vect>
    struct VectorCoordStruct<Vect, 0U> {
      
      using value_type = std::decay_t<decltype(std::declval<Vect>().X())>;
      
      static value_type get(const Vect& v) { return v.X(); }
      
      static Vect& set(Vect& v, value_type val) { v.SetX(val); return v;}
      
      static Vect& add(Vect& v, value_type val) { return set(v, get(v) + val); }
      
      static Vect& scale(Vect& v, value_type factor)
        { return set(v, get(v) * factor); }
      
    }; // VectorCoordStruct<0>
    
    template <typename Vect>
    struct VectorCoordStruct<Vect, 1U> {
      
      using value_type = std::decay_t<decltype(std::declval<Vect>().Y())>;
      
      static value_type get(const Vect& v) { return v.Y(); }
      
      static Vect& set(Vect& v, value_type val) { v.SetY(val); return v;}
      
      static Vect& add(Vect& v, value_type val) { return set(v, get(v) + val); }
      
      static Vect& scale(Vect& v, value_type factor)
        { return set(v, get(v) * factor); }
      
    }; // VectorCoordStruct<1>
    
    template <typename Vect>
    struct VectorCoordStruct<Vect, 2U> {
      
      using value_type = std::decay_t<decltype(std::declval<Vect>().Z())>;
      
      static value_type get(const Vect& v) { return v.Z(); }
      
      static Vect& set(Vect& v, value_type val) { v.SetZ(val); return v; }
      
      static Vect& add(Vect& v, value_type val) { return set(v, get(v) + val); }
      
      static Vect& scale(Vect& v, value_type factor)
        { return set(v, get(v) * factor); }
      
    }; // VectorCoordStruct<2>
    
    template <typename Vect>
    struct VectorCoordStruct<Vect, 3U> {
      
      using value_type = std::decay_t<decltype(std::declval<Vect>().T())>;
      
      static value_type get(const Vect& v) { return v.T(); }
      
      static Vect& set(Vect& v, value_type val) { v.SetT(val); return v; }
      
      static Vect& add(Vect& v, value_type val) { return set(v, get(v) + val); }
      
      static Vect& scale(Vect& v, value_type factor)
        { return set(v, get(v) * factor); }
      
    }; // VectorCoordStruct<3>
    
    
    template <typename Vect, unsigned int UpToDim = (VectorDim<Vect>::value - 1)>
    struct IsVectorHomogeneous
      : public std::integral_constant<bool,
          std::is_same<
            typename VectorCoordStruct<Vect, 0U>::value_type,
            typename VectorCoordStruct<Vect, UpToDim>::value_type
            >::value
          && IsVectorHomogeneous<Vect, (UpToDim - 1)>::value
        >
      {};
    
    template <typename Vect>
    struct IsVectorHomogeneous<Vect, 0U>: public std::true_type {};
    
    template <typename RefVect, typename... Vects>
    struct AreVectorsWithSameDimension;
    
    template <typename RefVect, typename First, typename... Others>
    struct AreVectorsWithSameDimension<RefVect, First, Others...>
      : public std::integral_constant<bool,
          (VectorDim<RefVect>::value == VectorDim<First>::value)
          && AreVectorsWithSameDimension<RefVect, Others...>::value
        >
      {};
    
    template <typename RefVect, typename First>
    struct AreVectorsWithSameDimension<RefVect, First>
      : public std::true_type {};
    
    template <typename Vect>
    using VectorValue_t = std::enable_if_t<
      IsVectorHomogeneous<Vect>::value,
      typename VectorCoordStruct<Vect, 0U>::value_type
      >;
    
    
    
    template <unsigned int Dim, typename Vect>
    auto VectorCoord(const Vect& v)
      { return VectorCoordStruct<Vect, Dim>::get(v); }
    
    template <unsigned int Dim, typename Vect, typename T>
    auto& SetVectorCoord(Vect& v, T val)
      { return VectorCoordStruct<Vect, Dim>::set(v, val); }
    
    template <unsigned int Dim, typename Vect, typename T>
    auto& AddToVectorCoord(Vect& v, T val)
      { return VectorCoordStruct<Vect, Dim>::add(v, val); }
    
    template <unsigned int Dim, typename Vect, typename T>
    auto& ScaleVectorCoord(Vect& v, T factor)
      { return VectorCoordStruct<Vect, Dim>::scale(v, factor); }
    
    
    template<
      typename DestVect, typename SrcVect,
      unsigned int Dim = VectorDim<DestVect>::value - 1
      >
    struct AddToVector {
      static DestVect& add(DestVect& dest, SrcVect const& src)
        {
          AddToVector<DestVect, SrcVect, Dim-1>::add(dest, src);
          return AddToVectorCoord<Dim>(dest, VectorCoord<Dim>(src));
        }
    }; // struct AddToVector
    
    template<typename DestVect, typename SrcVect>
    struct AddToVector<DestVect, SrcVect, 0U> {
      static DestVect& add(DestVect& dest, SrcVect const& src)
        { return AddToVectorCoord<0U>(dest, VectorCoord<0U>(src)); }
    }; // struct AddToVector<0>
    
    
    template
      <typename Vect, typename T, unsigned int Dim = (VectorDim<Vect>::value-1)>
    struct ScaleVectorCoords {
      static Vect& scale(Vect& v, T factor)
        {
          ScaleVectorCoords<Vect, T, Dim-1>::scale(v, factor);
          return ScaleVectorCoord<Dim>(v, factor);
        }
    }; // ScaleVectorCoords
    
    template <typename Vect, typename T>
    struct ScaleVectorCoords<Vect, T, 0U> {
      static Vect& scale(Vect& v, T factor)
        { return ScaleVectorCoord<0U>(v, factor); }
    }; // ScaleVectorCoords<0>
    
    
    template <typename DestVect, typename SrcVect>
    DestVect& addToVector(DestVect& dest, SrcVect const& src)
      { return AddToVector<DestVect, SrcVect>::add(dest, src); }
    
    template <typename Vect, typename T>
    Vect& scaleVector(Vect& v, T factor)
      { return ScaleVectorCoords<Vect, T>::scale(v, factor); }
    
    //--------------------------------------------------------------------------
    
    template <typename DestVect, typename... SrcVects>
    struct AccumulateVectors;
      
    template <typename DestVect, typename First, typename... Others>
    struct AccumulateVectors<DestVect, First, Others...> {
      static DestVect& add
        (DestVect& sum, First const& first, Others const&... others)
        {
          addToVector(sum, first);
          return AccumulateVectors<DestVect, Others...>::add(sum, others...);
        }
    };
    
    template <typename DestVect>
    struct AccumulateVectors<DestVect> {
      static DestVect& add(DestVect& sum) { return sum; }
    };
    
    
    template <typename... Vects>
    struct VectorSummer;
    
    template <typename First, typename... Others>
    struct VectorSummer<First, Others...> {
      static First sum(First const& first, Others const&... others)
        {
          using DestVect = First;
          DestVect s = first;
          return AccumulateVectors<DestVect, Others...>::add(s, others...);
        }
    }; // struct VectorSummer
    
    
    template <typename... Vects>
    auto sumVectors(Vects const&... vs)
      {
        static_assert
          (sizeof...(Vects) > 0, "sumVectors requires at least one argument");
        return VectorSummer<Vects...>::sum(vs...);
      } // sumVectors()
    
    
    namespace is_iterable_impl_ns {
      // we use a namespace to allow lookup of STL begin() and end()
      // while still allowing for custom versions (this is fairly standard)
      using std::begin;
      using std::end;
      
      // when asking for is_iterable_impl<T>(0), is_iterable_impl(int)
      // is preferred if available; it is available if the expression in
      // decltype() is valid. The expression uses comma operator, and the result
      // of decltype() expression is the one of the last branch of it, that is
      // true_type, regardless of the result of the other branches -- just the
      // usual comma operator behaviour; unfortunately comma operators may be
      // overloaded, so we interleave void() elements so that the operators are
      // operator,(bool, void) (which takes basic types and can't be overloaded)
      // just in case. Back to the expression: for it to be valid, begin(T) and
      // end(T) need to be somehow available and of types that can be compared,
      // and the comparison result converted to boolean.
      template <typename T>
      auto is_iterable_impl(int)
        -> decltype (
          bool(begin(std::declval<T&>()) != end(std::declval<T&>())),
          void(),
          std::true_type()
        );
      
      // this is the overload of is_iterable_impl() that gets picked whenever
      // the previous one is not available (or the argument is not integer)
      template <typename T>
      std::false_type is_iterable_impl(...);
      
    } // is_iterable_impl_ns
    
    template <typename T>
    using is_iterable = decltype(is_iterable_impl_ns::is_iterable_impl<T>(0));
    
    
    template <typename... Args>
    struct first_type;
    
    template <typename First, typename... Others>
    struct first_type<First, Others...> {
      using type = First;
    }; // first_type
    
    template <typename... Args>
    using first_type_t = typename first_type<Args...>::type;
    
    
    //--------------------------------------------------------------------------
    
    template <typename ... Vects>
    struct MiddlePointForVectors {
      static auto run(Vects const&... vs)
        {
          // determine the types involved
          using sum_type = decltype(sumVectors(vs...));
          using value_type = VectorValue_t<sum_type>;
          
          auto sum = sumVectors(vs...);
          return scaleVector(sum, value_type(1.) / sizeof...(vs));
        }
    }; // struct MiddlePointForVectors
    
    template <typename Coll>
    struct MiddlePointForCollection {
      static auto run(Coll const& vs)
        {
          // determine the types involved
          using sum_type = std::decay_t<decltype(*(vs.begin()))>;
          using value_type = VectorValue_t<sum_type>;
          
          assert(!vs.empty());
          
          using std::begin;
          using std::end;
          auto it = begin(vs);
          
          sum_type sum = *it;
          if (vs.size() == 1) return sum;
          
          auto vend = end(vs);
          while (++it != vend) addToVector(sum, *it);
          return scaleVector(sum, value_type(1.) / vs.size());
        }
    }; // struct MiddlePointForCollection
    
    
    //--------------------------------------------------------------------------
    
    template <
      typename Stream, typename Vect,
      unsigned int UpToDim = (VectorDim<Vect>::value - 1)
      >
    struct VectorPrinter {
      static_assert(VectorDim<Vect>::value >= 1,
        "Only vectors with dimension larger than 0 can be printed");
      static_assert(UpToDim < VectorDim<Vect>::value,
        "Can't print vector beyond its dimensions");
      static void print(Stream&& out, Vect const& v)
        {
          VectorPrinter<Stream, Vect, (UpToDim - 1)>::print
            (std::forward<Stream>(out), v);
          out << ", " << VectorCoord<UpToDim>(v);
        }
    }; // VectorPrinter
    
    template <typename Stream, typename Vect>
    struct VectorPrinter<Stream, Vect, 0U>
    {
      static void print(Stream&& out, Vect const& v)
        { out << VectorCoord<0U>(v); }
    }; // VectorPrinter<0>
    
    
    template <typename Stream, typename Vect>
    void printVectorCoords(Stream&& out, Vect const& v)
      {
        out << " ";
        if (VectorDim<Vect>::value < 1) return; // actually won't instantiate...
        VectorPrinter<Stream, Vect>::print(std::forward<Stream>(out), v);
        out << " ";
      }
  } // namespace details
  
  
  /// Returns the middle point among the specified vectors;
  /// the returned type is the one of the first argument.
  /// Each of the vectors has to support accessors X(), Y() and, for 3D, Z();
  /// the return type must also be constructible via its coordinates.
  template <typename First, typename... Others>
  auto middlePoint(First const& first, Others const&... others)
    {
      using MiddlePointAlg_t = std::conditional_t
        <
          ((sizeof...(Others) == 0) && details::is_iterable<First>::value)
          , details::MiddlePointForCollection<First>
          , details::MiddlePointForVectors<First, Others...>
        >;
      return MiddlePointAlg_t::run(first, others...);
    } // middlePoint()
  
  
  template <typename Vect>
  class PointAverage {
      public:
    using vector_type = Vect;
    using value_type = details::VectorValue_t<Vect>;
    
    
    void clear() { count = 0; }
    
    template <typename Other>
    void add(Other const& v)
      {
        if (count == 0) sum = v;
        else details::addToVector(sum, v);
        count += 1;
      }
    
    value_type n() const { return count; }
    
    Vect&& yield()
      {
        details::scaleVector(sum, value_type(1.) / count);
        count = 0;
        return std::move(sum);
      }
    
    
      private:
    value_type count = 0;
    vector_type sum;
    
  }; // class PointAverage
} // namespace pma


#endif // LARRECO_RECOALG_PMALG_GEOMDEFS_H
