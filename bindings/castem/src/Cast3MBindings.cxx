//
// Created by ds261764 on 27/01/2021.
//

#include <iostream>
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
#include <cstring>
#include <stdexcept>
#include "Mechhano/Cast3M/Cast3MBindings.hxx"
//#include "Eigen/Dense"
#include "lib/eigen3/Eigen/Dense"

extern "C" {

    #if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
      // code retrieved from
      // http://www.codeproject.com/Tips/479880/GetLastError-as-std-string
      static std::string getLastWin32Error() {
        const DWORD error = GetLastError();
        if (error) {
          LPVOID lpMsgBuf;
          DWORD bufLen = FormatMessage(
              FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                  FORMAT_MESSAGE_IGNORE_INSERTS,
              nullptr, error, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
              (LPTSTR)&lpMsgBuf, 0, nullptr);
          if (bufLen) {
            LPCSTR lpMsgStr = (LPTSTR)lpMsgBuf;
            std::string result(lpMsgStr, lpMsgStr + bufLen);
            LocalFree(lpMsgBuf);
            return result;
          }
        }
        return std::string();
      }
    #endif /*  (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */

      static std::string getErrorMessage() {
    #if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
        return getLastWin32Error();
    #else
        const auto e = ::dlerror();
        if (e != nullptr) {
          return std::string(e);
        }
        return "";
    #endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
      }    // end of  getErrorMessage


    castem_hho_status castem_hho_report_success() {
      return {CASTEM_HHO_SUCCESS, nullptr};
    }  // end of castem_hho_status castem_hho_report_success

    castem_hho_status castem_hho_report_failure(const char *const e) {
      static thread_local char msg[512];
      ::strncpy(msg, e, 511);
      msg[511] = '\0';
      return {CASTEM_HHO_FAILURE, msg};
    }  // end of castem_hho_status castem_hho_report_failure

    castem_hho_status castem_hho_handle_cxx_exception() {
      try {
        throw;
      } catch (std::exception &e) {
        return castem_hho_report_failure(e.what());
      } catch (...) {
        return castem_hho_report_failure(
            "castem_hho_handle_cxx_exception: "
            "unknown exception");
      }
    }  // end of castem_hho_handle_cxx_exception

    castem_hho_status castem_hho_get_element_functions(
        castem_hho_element_functions *const d,
//        const castem_hho_element_geometry *const gd,
        const char *const l,
        const char *const f) {
      using fct_ptr = castem_hho_status (*)(castem_hho_element_functions *const
//                                const castem_hho_element_geometry *const
                                );
      std::cout << l << std::endl;
      std::cout << f << std::endl;
    #if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
      const auto lib  = LoadLibrary(TEXT(l));
      if (lib == nullptr) {
        return castem_hho_report_failure(getErrorMessage().c_str());
      }
      const auto fct = reinterpret_cast<fct_ptr>(::GetProcAddress(lib, f));
      if(fct == nullptr){
        return castem_hho_report_failure(getErrorMessage().c_str());
      }
      return fct(d, gd);
    #else
      const auto lib = ::dlopen(l, RTLD_NOW);
      if (lib == nullptr) {
        return castem_hho_report_failure(getErrorMessage().c_str());
      }
      const auto fct = reinterpret_cast<fct_ptr>(dlsym(lib, f));
      if (fct == nullptr) {
        return castem_hho_report_failure(getErrorMessage().c_str());
      }
//      return fct(d, gd);
        return fct(d);
    #endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    }  // end of castem_hho_get_element_description

    castem_hho_status castem_hho_get_element_description(
            castem_hho_element_description *const d,
            const char *const l,
            const char *const f) {
        using fct_ptr = castem_hho_status (*)(castem_hho_element_description *const);
        std::cout << l << std::endl;
        std::cout << f << std::endl;
    #if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
        const auto lib  = LoadLibrary(TEXT(l));
          if (lib == nullptr) {
            return castem_hho_report_failure(getErrorMessage().c_str());
          }
          const auto fct = reinterpret_cast<fct_ptr>(::GetProcAddress(lib, f));
          if(fct == nullptr){
            return castem_hho_report_failure(getErrorMessage().c_str());
          }
          return fct(d, gd);
    #else
        const auto lib = ::dlopen(l, RTLD_NOW);
        if (lib == nullptr) {
            return castem_hho_report_failure(getErrorMessage().c_str());
        }
        const auto fct = reinterpret_cast<fct_ptr>(dlsym(lib, f));
        if (fct == nullptr) {
            return castem_hho_report_failure(getErrorMessage().c_str());
        }
        return fct(d);
    #endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    }  // end of castem_hho_get_element_description

    castem_hho_status castem_hho_get_generic_functions(
            castem_hho_generic_functions *const d,
            const char *const l,
            const char *const f
            ) {
        using fct_ptr = castem_hho_status (*)(castem_hho_generic_functions *const);
        std::cout << l << std::endl;
        std::cout << f << std::endl;
    #if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
        const auto lib  = LoadLibrary(TEXT(l));
              if (lib == nullptr) {
                return castem_hho_report_failure(getErrorMessage().c_str());
              }
              const auto fct = reinterpret_cast<fct_ptr>(::GetProcAddress(lib, f));
              if(fct == nullptr){
                return castem_hho_report_failure(getErrorMessage().c_str());
              }
              return fct(d, gd);
    #else
        const auto lib = ::dlopen(l, RTLD_NOW);
        if (lib == nullptr) {
            return castem_hho_report_failure(getErrorMessage().c_str());
        }
        const auto fct = reinterpret_cast<fct_ptr>(dlsym(lib, f));
        if (fct == nullptr) {
            return castem_hho_report_failure(getErrorMessage().c_str());
        }
        return fct(d);
    #endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    }  // end of castem_hho_get_element_description

    castem_hho_status castem_hho_get_gradient_operator(
            const castem_hho_element_functions *const elem_funs,
            const castem_hho_element_geometry *const elem_geom,
            double *const data, // worksapce
            int64_t index
    ) {
        return elem_funs->get_gradient_operator(
                elem_geom,
                data,
                index
        );
    }

    castem_hho_status castem_hho_get_stabilization_operator(
            const castem_hho_element_functions *const elem_funs,
            const castem_hho_element_geometry *const elem_geom,
            double *const data // worksapce
    ) {
        return elem_funs->get_stabilization_operator(
                elem_geom,
                data
        );
    }

    castem_hho_status castem_hho_get_gauss_weight(
            const castem_hho_element_functions *const elem_funs,
            const castem_hho_element_geometry *const elem_geom,
            double *const data, // worksapce
            int64_t index
    ) {
        return elem_funs->get_gauss_weight(
                elem_geom,
                data,
                index
        );
    }

    castem_hho_status castem_hho_get_gauss_point(
            const castem_hho_element_functions *const elem_funs,
            const castem_hho_element_geometry *const elem_geom,
            double *const data, // worksapce
            int64_t index
    ) {
        return elem_funs->get_gauss_point(
                elem_geom,
                data,
                index
        );
    }

    castem_hho_status castem_hho_invert_matrix(
            const castem_hho_generic_functions *const gene_funs,
            double *const data, // data
            int64_t index // index
            ) {
        return gene_funs->invert_matrix(
                data,
                index
        );
    }

}  // end of extern "C"
