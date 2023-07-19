#ifndef MACROS_HPP
#define MACROS_HPP

/*! \file macros.hpp
 *  \brief This file contains macros used throughout
 */

/// @brief computes the square of its argument
#define SQUARE(x) ((x)*(x))

/// @brief define a private variable and a public getter / setter pair
#define METHLASSO_GET_SET_FULL_DESCR(GType, SType, MType, Name, Descr)         \
public:                                                          \
  /*! @{ */                                                      \
  /** @brief Gets and sets Name                                 \
   Descr                                                         \
   */                                                            \
  GType get_##Name() const {                                      \
    return Name##_;                                              \
  };                                                             \
  void set_##Name(SType value) {                                 \
    Name##_ = value;                                             \
  }                                                              \
  /*! @} */                                                      \
  private:                                                       \
    MType Name##_;

  
/// @brief define a private variable and a public getter / setter pair, set/return by value
#define METHLASSO_GET_SET_SIMPLE_DESCR(Type, Name, Descr)         \
  METHLASSO_GET_SET_FULL_DESCR(Type, Type, Type, Name, Descr)

#define METHLASSO_GET_SET_SIMPLE(Type, Name)         \
  METHLASSO_GET_SET_SIMPLE_DESCR(Type, Name,)


/// @brief define a private variable and a public getter / setter pair, set by const ref
#define METHLASSO_GET_SET_CONSTREF_DESCR(Type, Name, Descr)         \
  METHLASSO_GET_SET_FULL_DESCR(Type, const Type&, Type, Name, Descr)

#define METHLASSO_GET_SET_CONSTREF(Type, Name)         \
  METHLASSO_GET_SET_CONSTREF_DESCR(Type, Name,)


/// @brief define a private variable and a public getter / setter pair, set and return by const ref
#define METHLASSO_GET_SET_OWN_DESCR(Type, Name, Descr)         \
  METHLASSO_GET_SET_FULL_DESCR(const Type&, const Type&, Type, Name, Descr)

#define METHLASSO_GET_SET_OWN(Type, Name)         \
  METHLASSO_GET_SET_OWN_DESCR(Type, Name,)


/// @brief define a private variable and a public getter
#define METHLASSO_GET_FULL_DESCR(GType, MType, Name, Descr)         \
  public:                                                          \
  /*! @{ */                                                      \
  /** @brief Gets Name by const reference              \
   Descr                                                         \
   */                                                            \
    GType get_##Name() const {                               \
      return Name##_;                                              \
    };                                                             \
  /*! @} */                                                      \
    private:                                                       \
      MType Name##_;

        
/// @brief define a private const variable and a public getter, return by const ref
#define METHLASSO_GET_CONST_OWN_DESCR(Type, Name, Descr)              \
  METHLASSO_GET_FULL_DESCR(const Type&, const Type, Name, Descr)
  
#define METHLASSO_GET_CONST_OWN(Type, Name)              \
  METHLASSO_GET_CONST_OWN_DESCR(Type, Name,)


/// @brief define a private const variable and a public getter, return by value
#define METHLASSO_GET_CONST_SIMPLE_DESCR(Type, Name, Descr)              \
METHLASSO_GET_FULL_DESCR(Type, const Type, Name, Descr)

#define METHLASSO_GET_CONST_SIMPLE(Type, Name)              \
METHLASSO_GET_CONST_SIMPLE_DESCR(Type, Name,)


#endif

