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

    castem_hho_status castem_hho_print_element_geometry(
        const castem_hho_element_geometry *const gd) {
      if (gd == nullptr) {
        return castem_hho_report_failure(
            "castem_hho_print_element_geometry: "
            "invalid argument");
      }
      std::cout << "- euclidian dimension: " << gd->dim_eucli << '\n';
      std::cout << "- number of vertices_coordinates: " << gd->num_vertices << '\n';
      for (int i = 0; i != gd->num_vertices; ++i) {
        std::cout << "\t- vertex " << i << ": ";
        for (int d = 0; d != gd->dim_eucli; ++d) {
          std::cout << " " << gd->vertices_coordinates[i + d * gd->num_vertices];
        }
        std::cout << '\n';
      }
      std::cout << "- number of faces: " << gd->num_faces << '\n';
      std::cout << "- number of vertices_coordinates per face:\n";
      for (int i = 0; i != gd->num_faces; ++i) {
        std::cout << "\t- face " << i << ": " << gd->num_vertices_per_face[i] << '\n';
      }
      std::cout << "- faces connectivity:\n";
      auto pos = int{};
      for (int i = 0; i != gd->num_faces; ++i) {
        std::cout << "\t- face " << i << ":";
        for (int v = 0; v != gd->num_vertices_per_face[i]; ++v, ++pos) {
          std::cout << " " << gd->connectivity[pos];
        }
        std::cout << '\n';
      }
      return castem_hho_report_success();
    }; // end of castem_hho_print_element_geometry

    castem_hho_status castem_hho_get_element_description(
        castem_hho_element_description *const d,
        const char *const l,
        const char *const f,
        const castem_hho_element_geometry *const gd) {
      using fct_ptr =
          castem_hho_status (*)(castem_hho_element_description *const,
                                const castem_hho_element_geometry *const);
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
      return fct(d, gd);
    #endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    }  // end of castem_hho_get_element_description

    castem_hho_status castem_hho_initialize_workspace(
            const castem_hho_element_description *const ed,
            double *const workspace, // worksapce
            const castem_hho_element_geometry *const gd
    ) {
        return ed->initialize_workspace(
                workspace,
                gd
        );
    }

    castem_hho_status castem_hho_get_gauss_data(
            const castem_hho_element_description *const ed,
            double *const pts,
            double *const wts,
            const castem_hho_element_geometry *const gd
    ) {
        return ed->get_gauss_data(
                pts,
                wts,
                gd
        );
    }

    castem_hho_status castem_hho_compute_gradients(
            const castem_hho_element_description *const ed,
            double *const workspace, // worksapce
            double *const gradient_data, // OUT
            const double *const faces_unknows, // IN
            const castem_hho_element_geometry *const gd
            ) {
        return ed->compute_gradients(
                workspace,
                gradient_data,
                faces_unknows,
                gd
                );
    }

    castem_hho_status castem_hho_compute_internal_forces(
            const castem_hho_element_description *const ed,
            double *const workspace, // worksapce
            double *const internal_forces_data, // OUT
            const double *const stress, // IN
            const double *const faces_unknows, // IN
            const castem_hho_element_geometry *const gd
    ) {
        return ed->compute_internal_forces(
                workspace,
                internal_forces_data,
                stress,
                faces_unknows,
                gd
        );
    }

    castem_hho_status castem_hho_compute_system(
            const castem_hho_element_description *const ed,
            double *const workspace, // worksapce
            double *const system_data, // OUT
            const double *const ktan, // IN
            const double *const res, // IN
            const castem_hho_element_geometry *const gd
    ) {
        return ed->compute_system(
                workspace,
                system_data,
                ktan,
                res,
                gd
        );
    }

    castem_hho_status castem_hho_decondensate(
            const castem_hho_element_description *const ed,
            double *const workspace, // worksapce
            const double *const faces_corrections, // IN
            const castem_hho_element_geometry *const gd
    ) {
        return ed->decondensate(
                workspace,
                faces_corrections,
                gd
        );
    }

}  // end of extern "C"
