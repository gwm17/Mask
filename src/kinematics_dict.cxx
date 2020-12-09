// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIkinematics_dict
#define R__NO_DEPRECATION

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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/Kinematics.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *NucData_Dictionary();
   static void NucData_TClassManip(TClass*);
   static void *new_NucData(void *p = 0);
   static void *newArray_NucData(Long_t size, void *p);
   static void delete_NucData(void *p);
   static void deleteArray_NucData(void *p);
   static void destruct_NucData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NucData*)
   {
      ::NucData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NucData));
      static ::ROOT::TGenericClassInfo 
         instance("NucData", "include/Kinematics.h", 13,
                  typeid(::NucData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NucData_Dictionary, isa_proxy, 4,
                  sizeof(::NucData) );
      instance.SetNew(&new_NucData);
      instance.SetNewArray(&newArray_NucData);
      instance.SetDelete(&delete_NucData);
      instance.SetDeleteArray(&deleteArray_NucData);
      instance.SetDestructor(&destruct_NucData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NucData*)
   {
      return GenerateInitInstanceLocal((::NucData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NucData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NucData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::NucData*)0x0)->GetClass();
      NucData_TClassManip(theClass);
   return theClass;
   }

   static void NucData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_NucData(void *p) {
      return  p ? new(p) ::NucData : new ::NucData;
   }
   static void *newArray_NucData(Long_t nElements, void *p) {
      return p ? new(p) ::NucData[nElements] : new ::NucData[nElements];
   }
   // Wrapper around operator delete
   static void delete_NucData(void *p) {
      delete ((::NucData*)p);
   }
   static void deleteArray_NucData(void *p) {
      delete [] ((::NucData*)p);
   }
   static void destruct_NucData(void *p) {
      typedef ::NucData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NucData

namespace {
  void TriggerDictionaryInitialization_kinematics_dict_Impl() {
    static const char* headers[] = {
"include/Kinematics.h",
0
    };
    static const char* includePaths[] = {
"/home/gordon/cern/root-6.22.02/root-install/include/",
"/home/gordon/Kinematics/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "kinematics_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate("$clingAutoload$include/Kinematics.h")))  NucData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "kinematics_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "include/Kinematics.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"NucData", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("kinematics_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_kinematics_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_kinematics_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_kinematics_dict() {
  TriggerDictionaryInitialization_kinematics_dict_Impl();
}
