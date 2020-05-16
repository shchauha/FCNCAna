// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LinkDef_out

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/QuantFuncMathCore.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/TLorentzVector.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/Vector4D.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR_Dictionary();
   static void ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR_TClassManip(TClass*);
   static void *new_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p = 0);
   static void *newArray_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(Long_t size, void *p);
   static void delete_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p);
   static void deleteArray_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p);
   static void destruct_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)
   {
      ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", "Math/GenVector/LorentzVector.h", 48,
                  typeid(::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >) );
      instance.SetNew(&new_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR);
      instance.SetNewArray(&newArray_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR);
      instance.SetDelete(&delete_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR);
      instance.SetDestructor(&destruct_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)
   {
      return GenerateInitInstanceLocal((::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)0x0)->GetClass();
      ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR_Dictionary();
   static void ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR_TClassManip(TClass*);
   static void *new_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p = 0);
   static void *newArray_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(Long_t size, void *p);
   static void delete_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p);
   static void deleteArray_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p);
   static void destruct_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT::Math::PxPyPzE4D<float>*)
   {
      ::ROOT::Math::PxPyPzE4D<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ROOT::Math::PxPyPzE4D<float>));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::Math::PxPyPzE4D<float>", "Math/GenVector/PxPyPzE4D.h", 42,
                  typeid(::ROOT::Math::PxPyPzE4D<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(::ROOT::Math::PxPyPzE4D<float>) );
      instance.SetNew(&new_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR);
      instance.SetNewArray(&newArray_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR);
      instance.SetDelete(&delete_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR);
      instance.SetDestructor(&destruct_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT::Math::PxPyPzE4D<float>*)
   {
      return GenerateInitInstanceLocal((::ROOT::Math::PxPyPzE4D<float>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ROOT::Math::PxPyPzE4D<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ROOT::Math::PxPyPzE4D<float>*)0x0)->GetClass();
      ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > : new ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >;
   }
   static void *newArray_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >[nElements] : new ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p) {
      delete ((::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)p);
   }
   static void deleteArray_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p) {
      delete [] ((::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >*)p);
   }
   static void destruct_ROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgR(void *p) {
      typedef ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::ROOT::Math::PxPyPzE4D<float> : new ::ROOT::Math::PxPyPzE4D<float>;
   }
   static void *newArray_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::ROOT::Math::PxPyPzE4D<float>[nElements] : new ::ROOT::Math::PxPyPzE4D<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p) {
      delete ((::ROOT::Math::PxPyPzE4D<float>*)p);
   }
   static void deleteArray_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p) {
      delete [] ((::ROOT::Math::PxPyPzE4D<float>*)p);
   }
   static void destruct_ROOTcLcLMathcLcLPxPyPzE4DlEfloatgR(void *p) {
      typedef ::ROOT::Math::PxPyPzE4D<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ROOT::Math::PxPyPzE4D<float>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 214,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR_Dictionary();
   static void vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p = 0);
   static void *newArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p);
   static void deleteArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p);
   static void destruct_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >*)
   {
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >", -2, "vector", 214,
                  typeid(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >) );
      instance.SetNew(&new_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >*)0x0)->GetClass();
      vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > : new vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >;
   }
   static void *newArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >[nElements] : new vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p) {
      delete ((vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >*)p);
   }
   static void deleteArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p) {
      delete [] ((vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >*)p);
   }
   static void destruct_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEfloatgRsPgRsPgR(void *p) {
      typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >

namespace {
  void TriggerDictionaryInitialization_LinkDef_out_Impl() {
    static const char* headers[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/QuantFuncMathCore.h",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/TLorentzVector.h",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/Vector4D.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2/include",
"/home/users/fgolf/fcnc/babymaking/batch/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LinkDef_out dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace ROOT{namespace Math{template <class ScalarType = double> class __attribute__((annotate("$clingAutoload$Math/GenVector/PxPyPzE4D.h")))  __attribute__((annotate("$clingAutoload$/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/Vector4D.h")))  PxPyPzE4D;
}}
namespace ROOT{namespace Math{template <class CoordSystem> class __attribute__((annotate("$clingAutoload$Math/Vector4Dfwd.h")))  __attribute__((annotate("$clingAutoload$/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/Vector4D.h")))  LorentzVector;
}}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LinkDef_out dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/QuantFuncMathCore.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/TLorentzVector.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike2//include/Math/Vector4D.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", payloadCode, "@",
"ROOT::Math::PxPyPzE4D<float>", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LinkDef_out",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LinkDef_out_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LinkDef_out_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LinkDef_out() {
  TriggerDictionaryInitialization_LinkDef_out_Impl();
}
