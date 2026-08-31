// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME simEventDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "G4d2oGeom.h"
#include "simEvent.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void delete_G4d2oGeom(void *p);
   static void deleteArray_G4d2oGeom(void *p);
   static void destruct_G4d2oGeom(void *p);
   static void streamer_G4d2oGeom(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::G4d2oGeom*)
   {
      ::G4d2oGeom *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::G4d2oGeom >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("G4d2oGeom", ::G4d2oGeom::Class_Version(), "G4d2oGeom.h", 8,
                  typeid(::G4d2oGeom), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::G4d2oGeom::Dictionary, isa_proxy, 16,
                  sizeof(::G4d2oGeom) );
      instance.SetDelete(&delete_G4d2oGeom);
      instance.SetDeleteArray(&deleteArray_G4d2oGeom);
      instance.SetDestructor(&destruct_G4d2oGeom);
      instance.SetStreamerFunc(&streamer_G4d2oGeom);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::G4d2oGeom*)
   {
      return GenerateInitInstanceLocal(static_cast<::G4d2oGeom*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::G4d2oGeom*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_simHit(void *p = nullptr);
   static void *newArray_simHit(Long_t size, void *p);
   static void delete_simHit(void *p);
   static void deleteArray_simHit(void *p);
   static void destruct_simHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::simHit*)
   {
      ::simHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::simHit >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("simHit", ::simHit::Class_Version(), "simEvent.h", 7,
                  typeid(::simHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::simHit::Dictionary, isa_proxy, 4,
                  sizeof(::simHit) );
      instance.SetNew(&new_simHit);
      instance.SetNewArray(&newArray_simHit);
      instance.SetDelete(&delete_simHit);
      instance.SetDeleteArray(&deleteArray_simHit);
      instance.SetDestructor(&destruct_simHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::simHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::simHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::simHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_simAreaHit(void *p = nullptr);
   static void *newArray_simAreaHit(Long_t size, void *p);
   static void delete_simAreaHit(void *p);
   static void deleteArray_simAreaHit(void *p);
   static void destruct_simAreaHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::simAreaHit*)
   {
      ::simAreaHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::simAreaHit >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("simAreaHit", ::simAreaHit::Class_Version(), "simEvent.h", 28,
                  typeid(::simAreaHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::simAreaHit::Dictionary, isa_proxy, 4,
                  sizeof(::simAreaHit) );
      instance.SetNew(&new_simAreaHit);
      instance.SetNewArray(&newArray_simAreaHit);
      instance.SetDelete(&delete_simAreaHit);
      instance.SetDeleteArray(&deleteArray_simAreaHit);
      instance.SetDestructor(&destruct_simAreaHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::simAreaHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::simAreaHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::simAreaHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_simEvent(void *p = nullptr);
   static void *newArray_simEvent(Long_t size, void *p);
   static void delete_simEvent(void *p);
   static void deleteArray_simEvent(void *p);
   static void destruct_simEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::simEvent*)
   {
      ::simEvent *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::simEvent >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("simEvent", ::simEvent::Class_Version(), "simEvent.h", 53,
                  typeid(::simEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::simEvent::Dictionary, isa_proxy, 4,
                  sizeof(::simEvent) );
      instance.SetNew(&new_simEvent);
      instance.SetNewArray(&newArray_simEvent);
      instance.SetDelete(&delete_simEvent);
      instance.SetDeleteArray(&deleteArray_simEvent);
      instance.SetDestructor(&destruct_simEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::simEvent*)
   {
      return GenerateInitInstanceLocal(static_cast<::simEvent*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::simEvent*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr G4d2oGeom::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *G4d2oGeom::Class_Name()
{
   return "G4d2oGeom";
}

//______________________________________________________________________________
const char *G4d2oGeom::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::G4d2oGeom*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int G4d2oGeom::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::G4d2oGeom*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *G4d2oGeom::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::G4d2oGeom*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *G4d2oGeom::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::G4d2oGeom*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr simHit::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *simHit::Class_Name()
{
   return "simHit";
}

//______________________________________________________________________________
const char *simHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::simHit*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int simHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::simHit*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *simHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::simHit*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *simHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::simHit*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr simAreaHit::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *simAreaHit::Class_Name()
{
   return "simAreaHit";
}

//______________________________________________________________________________
const char *simAreaHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::simAreaHit*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int simAreaHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::simAreaHit*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *simAreaHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::simAreaHit*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *simAreaHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::simAreaHit*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr simEvent::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *simEvent::Class_Name()
{
   return "simEvent";
}

//______________________________________________________________________________
const char *simEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::simEvent*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int simEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::simEvent*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *simEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::simEvent*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *simEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::simEvent*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void G4d2oGeom::Streamer(TBuffer &R__b)
{
   // Stream an object of class G4d2oGeom.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_G4d2oGeom(void *p) {
      delete (static_cast<::G4d2oGeom*>(p));
   }
   static void deleteArray_G4d2oGeom(void *p) {
      delete [] (static_cast<::G4d2oGeom*>(p));
   }
   static void destruct_G4d2oGeom(void *p) {
      typedef ::G4d2oGeom current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_G4d2oGeom(TBuffer &buf, void *obj) {
      ((::G4d2oGeom*)obj)->::G4d2oGeom::Streamer(buf);
   }
} // end of namespace ROOT for class ::G4d2oGeom

//______________________________________________________________________________
void simHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class simHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(simHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(simHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_simHit(void *p) {
      return  p ? new(p) ::simHit : new ::simHit;
   }
   static void *newArray_simHit(Long_t nElements, void *p) {
      return p ? new(p) ::simHit[nElements] : new ::simHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_simHit(void *p) {
      delete (static_cast<::simHit*>(p));
   }
   static void deleteArray_simHit(void *p) {
      delete [] (static_cast<::simHit*>(p));
   }
   static void destruct_simHit(void *p) {
      typedef ::simHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::simHit

//______________________________________________________________________________
void simAreaHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class simAreaHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(simAreaHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(simAreaHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_simAreaHit(void *p) {
      return  p ? new(p) ::simAreaHit : new ::simAreaHit;
   }
   static void *newArray_simAreaHit(Long_t nElements, void *p) {
      return p ? new(p) ::simAreaHit[nElements] : new ::simAreaHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_simAreaHit(void *p) {
      delete (static_cast<::simAreaHit*>(p));
   }
   static void deleteArray_simAreaHit(void *p) {
      delete [] (static_cast<::simAreaHit*>(p));
   }
   static void destruct_simAreaHit(void *p) {
      typedef ::simAreaHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::simAreaHit

//______________________________________________________________________________
void simEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class simEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(simEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(simEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_simEvent(void *p) {
      return  p ? new(p) ::simEvent : new ::simEvent;
   }
   static void *newArray_simEvent(Long_t nElements, void *p) {
      return p ? new(p) ::simEvent[nElements] : new ::simEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_simEvent(void *p) {
      delete (static_cast<::simEvent*>(p));
   }
   static void deleteArray_simEvent(void *p) {
      delete [] (static_cast<::simEvent*>(p));
   }
   static void destruct_simEvent(void *p) {
      typedef ::simEvent current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::simEvent

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_libsimEvent() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_libsimEvent) = RecordReadRules_libsimEvent();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_libsimEvent));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_libsimEvent_Impl() {
    static const char* headers[] = {
"G4d2oGeom.h",
"simEvent.h",
nullptr
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/home/manoja450/geant4-install/include/Geant4",
"/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o",
"/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o/simEvent",
"/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o/include",
"/src",
"/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o/simEvent/include",
"/usr/include/root",
"/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o/include/simEvent/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libsimEvent dictionary forward declarations' payload"

#pragma diagnostic push
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class  __attribute__((annotate("$clingAutoload$G4d2oGeom.h")))  G4d2oGeom;
class  __attribute__((annotate("$clingAutoload$simEvent.h")))  simHit;
class  __attribute__((annotate("$clingAutoload$simEvent.h")))  simAreaHit;
class  __attribute__((annotate("$clingAutoload$simEvent.h")))  simEvent;
#pragma diagnostic pop
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libsimEvent dictionary payload"

#ifndef R__DUMMY_CXX_STANDARD_17
  #define R__DUMMY_CXX_STANDARD_17 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "G4d2oGeom.h"
#include "simEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"G4d2oGeom", payloadCode, "@",
"simAreaHit", payloadCode, "@",
"simEvent", payloadCode, "@",
"simHit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libsimEvent",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libsimEvent_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libsimEvent_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libsimEvent() {
  TriggerDictionaryInitialization_libsimEvent_Impl();
}
