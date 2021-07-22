//
// Created by dsiedel on 30/06/2021.
//

#ifndef MECHHCANO_FUNC_HXX
#define MECHHCANO_FUNC_HXX
#include <iostream>
#include "Cast3M/Cast3MBindings.hxx"

template<int a>
castem_hho_status initialize_workspace(
        double *const workspace, // worksapce
        const castem_hho_element_geometry *const gd
) {
    try {
        std::cout << "INITIALIZING WORKSPACE" << std::endl;
        workspace[0] = 3.0;
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

template<int a>
castem_hho_status get_gauss_data(
        double *const pts, // gauss points
        double *const wts, // gauss weights
        const castem_hho_element_geometry *const gd
) {
    try {
        std::cout << "WRITING GAUSS DATA" << std::endl;
        pts[0] = 13.6;
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

template<int a>
castem_hho_status compute_gradients(
        double *const workspace, //worksapce
        double *const gradient_data, // OUT
        const double *const faces_unknows, // IN
        const castem_hho_element_geometry *const gd
        ) {
    try {
        std::cout << "COMPUTING GRADIENTS" << std::endl;
        std::cout << "WK0 : " << workspace[0] << std::endl;
        gradient_data[0] = 1.0;
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

template<int a>
castem_hho_status compute_internal_forces(
        double *const workspace, //worksapce
        double *const internal_forces_data, // OUT
        const double *const stress, // IN
        const double *const faces_unknows, // IN
        const castem_hho_element_geometry *const gd
) {
    try {
        std::cout << "COMPUTING INTERNAL FORCES" << std::endl;
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

template<int a>
castem_hho_status compute_system(
        double *const workspace, //worksapce
        double *const system_data, // OUT
        const double *const ktan, // IN
        const double *const res, // IN
        const castem_hho_element_geometry *const gd
) {
    try {
        std::cout << "COMPUTING SYSTEM" << std::endl;
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

template<int a>
castem_hho_status getDescription(
        castem_hho_element_description *const d,
        const castem_hho_element_geometry *const gd
        ){
    try {
        std::cout << "INITIALIZING ELEMENT DESCRIPTION" << std::endl;
        d->element_workspace_size = 12;
        d->get_gauss_data = get_gauss_data<a>;
        d->initialize_workspace = initialize_workspace<a>;
        d->compute_gradients = compute_gradients<a>;
        d->compute_internal_forces = compute_internal_forces<a>;
        d->compute_system = compute_system<a>;
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

#endif //MECHHCANO_FUNC_HXX
